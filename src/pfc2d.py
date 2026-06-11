"""2D Phase-Field Crystal (Elder-Grant) solver, spectral semi-implicit scheme.

Free energy:  F = ∫ dx [ ψ/2 (r + (1+∇²)²) ψ + ψ⁴/4 ]
Dynamics:     ∂ψ/∂t = ∇² δF/δψ   (conserved)

Semi-implicit spectral update:
    ψ̂(t+dt) = [ψ̂(t) - dt k² FFT(ψ³)] / [1 + dt k² (r + (1-k²)²)]

Conventions follow BJK_hBN_Multiscale XPFC code (pyFFTW threads, periodic BC).
Tensile loading is applied by affine rescaling of the simulation box
(k-vector rescaling), following Hirouchi et al. / Stefanovic-style PFC
deformation: stretch x by (1+exx), contract y to conserve area by default.
"""

import os
import numpy as np

try:
    import pyfftw

    _NTHREADS = int(os.environ.get("PFC_FFT_THREADS", os.cpu_count() or 4))

    def _rfft2(a):
        return pyfftw.interfaces.numpy_fft.rfft2(a, threads=_NTHREADS)

    def _irfft2(a, s):
        return pyfftw.interfaces.numpy_fft.irfft2(a, s=s, threads=_NTHREADS)

    pyfftw.interfaces.cache.enable()
    FFT_BACKEND = "pyfftw"
except ImportError:
    _rfft2 = np.fft.rfft2

    def _irfft2(a, s):
        return np.fft.irfft2(a, s=s)

    FFT_BACKEND = "numpy"

# one-mode triangular lattice: q0 = 1, atomic spacing a0 = 4π/√3
A_LATTICE = 4.0 * np.pi / np.sqrt(3.0)


class PFC2D:
    def __init__(self, nx, ny, dx=np.pi / 4, r=-0.25, psi_bar=-0.25):
        self.nx, self.ny = nx, ny
        self.dx0 = dx          # undeformed grid spacing
        self.r = r
        self.psi_bar = psi_bar
        self.exx = 0.0         # accumulated true strain factors (box stretch)
        self.eyy = 0.0
        self.psi = np.full((ny, nx), psi_bar)
        self.time = 0.0
        self._update_k()

    # ---------- geometry / strain ----------
    @property
    def dx(self):
        return self.dx0 * (1.0 + self.exx)

    @property
    def dy(self):
        return self.dx0 * (1.0 + self.eyy)

    @property
    def lx(self):
        return self.nx * self.dx

    @property
    def ly(self):
        return self.ny * self.dy

    def _update_k(self):
        # real-FFT layout: last axis halved; Hermitian symmetry is then
        # structural, so no anti-Hermitian residue can grow in the k=1
        # amplification band (a complex-FFT version blew up after ~3000
        # uninterrupted steps from exactly that mode)
        kx = 2.0 * np.pi * np.fft.rfftfreq(self.nx, d=self.dx)
        ky = 2.0 * np.pi * np.fft.fftfreq(self.ny, d=self.dy)
        KX, KY = np.meshgrid(kx, ky)
        self.k2 = KX ** 2 + KY ** 2
        self.lin = self.r + (1.0 - self.k2) ** 2  # linear operator L(k)

    def apply_strain(self, dexx, area_conserving=True):
        """Increment box strain: stretch x by (1+dexx); optionally contract y
        so that (1+exx)(1+eyy) stays constant (area-conserving tension)."""
        new_fx = (1.0 + self.exx) * (1.0 + dexx)
        self.exx = new_fx - 1.0
        if area_conserving:
            self.eyy = 1.0 / (1.0 + self.exx) - 1.0
        self._update_k()

    # ---------- initial conditions ----------
    def init_random(self, noise=0.05, seed=0):
        rng = np.random.default_rng(seed)
        self.psi = self.psi_bar + noise * rng.standard_normal((self.ny, self.nx))

    def one_mode_field(self, X, Y, amp=None, ux=None, uy=None):
        """One-mode triangular solution with |k| = 1 wave vectors:
        ψ = ψ̄ + A [cos(qx x)cos(qy y) - 0.5 cos(2 qy y)], qx=√3/2, qy=1/2.
        qx, qy are snapped to the nearest box-commensurate values so the
        crystal is periodic across the boundaries.
        ux, uy: optional displacement fields (for seeding dislocations)."""
        r, pb = self.r, self.psi_bar
        if amp is None:
            # one-mode amplitude minimizing F (Elder & Grant 2004).
            # For psi_bar < 0 the triangular phase has A < 0: density maxima
            # then sit on the minima of f, which form a triangular lattice
            # (with A > 0 the maxima of f form a honeycomb instead).
            amp = -abs(pb + np.sqrt(-15.0 * r - 36.0 * pb ** 2) / 3.0) * 0.8
        qx, qy = self._snapped_q()
        Xe = X - (ux if ux is not None else 0.0)
        Ye = Y - (uy if uy is not None else 0.0)
        f = np.cos(qx * Xe) * np.cos(qy * Ye) - 0.5 * np.cos(2.0 * qy * Ye)
        return pb + amp * f * 4.0 / 3.0  # normalization of the 3-wave sum

    def init_crystal(self, noise=0.0, seed=0):
        x = np.arange(self.nx) * self.dx
        y = np.arange(self.ny) * self.dy
        X, Y = np.meshgrid(x, y)
        self.psi = self.one_mode_field(X, Y)
        if noise > 0:
            rng = np.random.default_rng(seed)
            self.psi += noise * rng.standard_normal((self.ny, self.nx))
        self._fix_mean()

    def _snapped_q(self):
        """Box-commensurate one-mode wave numbers (qx ~ √3/2, qy ~ 1/2)."""
        qx = 2.0 * np.pi * max(1, round(np.sqrt(3.0) / 2.0 * self.lx
                                        / (2 * np.pi))) / self.lx
        qy = 2.0 * np.pi * max(1, round(0.5 * self.ly / (2 * np.pi))) / self.ly
        return qx, qy

    def init_dislocations(self, cores, images=1):
        """Seed edge dislocations by phase winding (Skaugen et al. style):
        ux = (b/2π) Σ_images Σ_i sign_i · atan2(y - yi, x - xi), uy = 0.
        Each branch-cut jump equals one full Burgers vector (invisible to
        the lattice). `cores` is a list of (x_frac, y_frac, sign); signs
        must sum to zero for a periodic-compatible far field.

        cores=[(0.5, 0.25, +1), (0.5, 0.75, -1)] reproduces the classic
        vertical dipole; a quadrupole arrangement minimizes image
        interactions for multi-dislocation seeding."""
        if sum(s for _, _, s in cores) != 0:
            raise ValueError("net Burgers vector must be zero (signs sum 0)")
        x = np.arange(self.nx) * self.dx
        y = np.arange(self.ny) * self.dy
        X, Y = np.meshgrid(x, y)
        qx, _ = self._snapped_q()
        b = 2.0 * np.pi / qx  # x-period of the snapped lattice
        ux = np.zeros_like(X)
        for sx in range(-images, images + 1):
            for sy in range(-images, images + 1):
                for fx, fy, s in cores:
                    ux += s * np.arctan2(Y - fy * self.ly + sy * self.ly,
                                         X - fx * self.lx + sx * self.lx)
        ux *= b / (2.0 * np.pi)
        self.psi = self.one_mode_field(X, Y, ux=ux)
        self._fix_mean()

    def init_dislocation_dipole(self, images=1):
        """Classic +b/-b vertical dipole at (lx/2, ly/4) and (lx/2, 3ly/4)."""
        self.init_dislocations([(0.5, 0.25, +1), (0.5, 0.75, -1)],
                               images=images)

    def add_void(self, cx_frac=0.5, cy_frac=0.5, radius=None, edge=None,
                 depth=0.0):
        """Melt a disc as a stress concentrator (cf. the arc notch in the
        单晶铜拉伸模拟 project). With depth == 0 the disc is set to psi_bar
        and simply recrystallizes during relaxation (verified: useless as a
        notch). depth > 0 also depletes the local mean density
        (psi -> psi_bar - depth), putting the disc into the stable-liquid
        region; conserved dynamics then keeps the pore open like a real
        pore. Mass is intentionally NOT re-normalized."""
        if radius is None:
            radius = 3.0 * A_LATTICE
        if edge is None:
            edge = 0.5 * A_LATTICE
        x = np.arange(self.nx) * self.dx
        y = np.arange(self.ny) * self.dy
        X, Y = np.meshgrid(x, y)
        rr = np.sqrt((X - cx_frac * self.lx) ** 2 + (Y - cy_frac * self.ly) ** 2)
        w = 0.5 * (1.0 + np.tanh((rr - radius) / edge))  # 0 inside, 1 outside
        self.psi = w * self.psi + (1.0 - w) * (self.psi_bar - depth)
        if depth == 0.0:
            self._fix_mean()

    def _fix_mean(self):
        self.psi += self.psi_bar - self.psi.mean()

    # ---------- evolution ----------
    def step(self, dt, n=1, noise_amp=0.0, seed=None):
        """noise_amp > 0 adds approximately-conserved Gaussian fluctuations
        (zero-mean per step) — a thermal kick for nucleation studies, not a
        rigorous Model-B noise discretization."""
        rng = np.random.default_rng(seed) if noise_amp > 0 else None
        shape = self.psi.shape
        psi = self.psi
        # The real field is the state: a full rfft/irfft round trip EVERY
        # iteration projects out phantom spectral components. Carrying psi_h
        # across iterations lets anti-Hermitian residue (complex FFT: all
        # modes; rfft: the kx=0/Nyquist columns) grow unsaturated in the
        # k~1 amplification band -> box-spanning stripe blowup at fixed
        # physical time (~t=400 at 512^2). Verified empirically.
        for _ in range(n):
            if rng is not None:
                eta = rng.standard_normal(shape)
                psi = psi + np.sqrt(dt) * noise_amp * (eta - eta.mean())
            psi_h = _rfft2(psi)
            # stabilized splitting: C >= max 3 psi^2 keeps the explicit
            # remainder (psi^3 - C psi) contractive
            C = 3.0 * float(np.max(psi * psi))
            nl_h = _rfft2(psi ** 3)
            psi_h = ((psi_h - dt * self.k2 * (nl_h - C * psi_h))
                     / (1.0 + dt * self.k2 * (self.lin + C)))
            psi = _irfft2(psi_h, shape)
            self.time += dt
        self.psi = psi

    def stress(self, deps=1e-4):
        """Work-conjugate driving stress along the area-conserving tension
        path: sigma = d f / d exx at fixed psi (virtual affine deformation
        via k-rescaling; f = mean free-energy density)."""
        exx0, eyy0 = self.exx, self.eyy
        try:
            self.exx = (1.0 + exx0) * (1.0 + deps) - 1.0
            self.eyy = 1.0 / (1.0 + self.exx) - 1.0
            self._update_k()
            fp = self.free_energy()
            self.exx = (1.0 + exx0) * (1.0 - deps) - 1.0
            self.eyy = 1.0 / (1.0 + self.exx) - 1.0
            self._update_k()
            fm = self.free_energy()
        finally:
            self.exx, self.eyy = exx0, eyy0
            self._update_k()
        return (fp - fm) / (2.0 * deps * (1.0 + exx0))

    def free_energy(self):
        psi_h = _rfft2(self.psi)
        lin_term = _irfft2(self.lin * psi_h, self.psi.shape)
        f = 0.5 * self.psi * lin_term + 0.25 * self.psi ** 4
        return float(f.mean())

    # ---------- io ----------
    def save(self, path):
        np.savez_compressed(
            path, psi=self.psi, r=self.r, psi_bar=self.psi_bar,
            dx0=self.dx0, exx=self.exx, eyy=self.eyy, time=self.time)

    @classmethod
    def load(cls, path):
        d = np.load(path)
        m = cls(d["psi"].shape[1], d["psi"].shape[0], dx=float(d["dx0"]),
                r=float(d["r"]), psi_bar=float(d["psi_bar"]))
        m.psi = d["psi"]
        m.exx = float(d["exx"])
        m.eyy = float(d["eyy"])
        m.time = float(d["time"])
        m._update_k()
        return m
