"""3D Phase-Field Crystal solver (one-mode -> BCC), spectral semi-implicit.

Same free energy/dynamics as pfc2d, with 3D Laplacian. The one-mode BCC
solution uses the six {110} wave vectors:
    psi = psi_bar + A [cos(qx)cos(qy) + cos(qy)cos(qz) + cos(qz)cos(qx)],
q = 1/sqrt(2) so |k| = 1; BCC lattice constant a = 2*pi*sqrt(2).
Uniaxial tension: volume-conserving box rescaling along x
(eyy = ezz = 1/sqrt(1+exx) - 1).

FFT thread count respects env PFC_FFT_THREADS (to coexist with background
runs); defaults to all cores.
"""

import os
import numpy as np

try:
    import pyfftw

    _NTHREADS = int(os.environ.get("PFC_FFT_THREADS", os.cpu_count() or 4))

    def _rfftn(a):
        return pyfftw.interfaces.numpy_fft.rfftn(a, threads=_NTHREADS)

    def _irfftn(a, s):
        return pyfftw.interfaces.numpy_fft.irfftn(a, s=s, threads=_NTHREADS)

    pyfftw.interfaces.cache.enable()
    FFT_BACKEND = "pyfftw"
except ImportError:
    try:
        from scipy import fft as _sfft

        _WORKERS = int(os.environ.get("PFC_FFT_THREADS", os.cpu_count() or 4))

        def _rfftn(a):
            return _sfft.rfftn(a, workers=_WORKERS)

        def _irfftn(a, s):
            return _sfft.irfftn(a, s=s, workers=_WORKERS)

        FFT_BACKEND = "scipy"
    except ImportError:
        _rfftn = np.fft.rfftn

        def _irfftn(a, s):
            return np.fft.irfftn(a, s=s)

        FFT_BACKEND = "numpy"

A_BCC = 2.0 * np.pi * np.sqrt(2.0)   # BCC lattice constant for q0 = 1


class PFC3D:
    def __init__(self, nx, ny, nz, dx=np.pi / 4, r=-0.25, psi_bar=-0.25):
        self.nx, self.ny, self.nz = nx, ny, nz
        self.dx0 = dx
        self.r = r
        self.psi_bar = psi_bar
        self.exx = self.eyy = self.ezz = 0.0
        self.gxz = 0.0         # simple shear x += gxz*z (resolved across z-GB)
        self.psi = np.full((nz, ny, nx), psi_bar)
        self.time = 0.0
        self._update_k()

    @property
    def dx(self):
        return self.dx0 * (1.0 + self.exx)

    @property
    def dy(self):
        return self.dx0 * (1.0 + self.eyy)

    @property
    def dz(self):
        return self.dx0 * (1.0 + self.ezz)

    @property
    def lx(self):
        return self.nx * self.dx

    @property
    def ly(self):
        return self.ny * self.dy

    @property
    def lz(self):
        return self.nz * self.dz

    def _update_k(self):
        # real-FFT layout (last axis halved) — see pfc2d._update_k
        kx = 2 * np.pi * np.fft.rfftfreq(self.nx, d=self.dx)
        ky = 2 * np.pi * np.fft.fftfreq(self.ny, d=self.dy)
        kz = 2 * np.pi * np.fft.fftfreq(self.nz, d=self.dz)
        KZ, KY, KX = np.meshgrid(kz, ky, kx, indexing="ij")
        # simple xz shear F=[[1,0,g],[0,1,0],[0,0,1]] -> F^-T k = (kx,ky,kz-g*kx)
        KZP = KZ - self.gxz * KX
        self.k2 = KX ** 2 + KY ** 2 + KZP ** 2
        self.lin = self.r + (1.0 - self.k2) ** 2

    def apply_strain(self, dexx, volume_conserving=True):
        new_fx = (1.0 + self.exx) * (1.0 + dexx)
        self.exx = new_fx - 1.0
        if volume_conserving:
            f_t = 1.0 / np.sqrt(1.0 + self.exx)
            self.eyy = self.ezz = f_t - 1.0
        self._update_k()

    def apply_shear(self, dg):
        """Simple xz shear gamma_xz += dg (x += gamma*z). Applies resolved
        shear across a (001) twist GB; volume-preserving and free of the
        Bain/amorphization path of uniaxial BCC tension (the ~10-12% ceiling
        from M17), so it can probe GB glide/emission to larger strain."""
        self.gxz += dg
        self._update_k()

    def shear_stress(self, dg=1e-4):
        g0 = self.gxz
        try:
            self.gxz = g0 + dg
            self._update_k()
            fp = self.free_energy()
            self.gxz = g0 - dg
            self._update_k()
            fm = self.free_energy()
        finally:
            self.gxz = g0
            self._update_k()
        return (fp - fm) / (2.0 * dg)

    # ---------- initial conditions ----------
    def init_random(self, noise=0.05, seed=0):
        rng = np.random.default_rng(seed)
        self.psi = self.psi_bar + noise * rng.standard_normal(
            (self.nz, self.ny, self.nx))

    def _snapped_q(self):
        q = 1.0 / np.sqrt(2.0)
        qx = 2 * np.pi * max(1, round(q * self.lx / (2 * np.pi))) / self.lx
        qy = 2 * np.pi * max(1, round(q * self.ly / (2 * np.pi))) / self.ly
        qz = 2 * np.pi * max(1, round(q * self.lz / (2 * np.pi))) / self.lz
        return qx, qy, qz

    def init_crystal(self, amp=None):
        """One-mode BCC. Unlike the 2D triangular case, the BCC branch has
        A > 0 for psi_bar < 0: density maxima sit on BCC sites (0,0,0) and
        (pi,pi,pi) where f = 3 (verified energetically: A=+0.25 relaxes to
        F=0.0212 vs A=-0.25 at 0.0244 with 3x spurious peaks)."""
        if amp is None:
            amp = 0.25
        x = np.arange(self.nx) * self.dx
        y = np.arange(self.ny) * self.dy
        z = np.arange(self.nz) * self.dz
        Z, Y, X = np.meshgrid(z, y, x, indexing="ij")
        qx, qy, qz = self._snapped_q()
        f = (np.cos(qx * X) * np.cos(qy * Y)
             + np.cos(qy * Y) * np.cos(qz * Z)
             + np.cos(qz * Z) * np.cos(qx * X))
        self.psi = self.psi_bar + amp * f
        self._fix_mean()

    def _bcc_field(self, X, Y, Z, amp=0.25):
        qx, qy, qz = self._snapped_q()
        f = (np.cos(qx * X) * np.cos(qy * Y)
             + np.cos(qy * Y) * np.cos(qz * Z)
             + np.cos(qz * Z) * np.cos(qx * X))
        return self.psi_bar + amp * f

    def init_dislocation_lines(self, lines, amp=0.25):
        """Seed straight dislocation lines by 3D phase winding. Each line:
        dict(axis='x'|'y'|'z', pos=(a,b) [the two transverse fractional coords],
             burgers='x'|'y'|'z', sign=+1|-1). The winding angle theta rotates
        in the plane perpendicular to the line axis about (pos); the
        displacement is along the Burgers direction, u = (b/2pi)*sign*theta.
        Net Burgers per component must cancel (PBC). Forest = lines on
        different axes so they thread each other's glide planes -> junction
        tests (the 3D process 2D PFC cannot do)."""
        x = np.arange(self.nx) * self.dx
        y = np.arange(self.ny) * self.dy
        z = np.arange(self.nz) * self.dz
        Z, Y, X = np.meshgrid(z, y, x, indexing="ij")
        qx, _, _ = self._snapped_q()
        bmag = 2.0 * np.pi / qx
        coords = {"x": (X, self.lx), "y": (Y, self.ly), "z": (Z, self.lz)}
        U = {"x": np.zeros_like(X), "y": np.zeros_like(X), "z": np.zeros_like(X)}
        net = {"x": 0.0, "y": 0.0, "z": 0.0}
        for ln in lines:
            ax = ln["axis"]
            transv = [a for a in ("x", "y", "z") if a != ax]
            (C0, L0), (C1, L1) = coords[transv[0]], coords[transv[1]]
            p0, p1 = ln["pos"]
            theta = np.arctan2(C1 - p1 * L1, C0 - p0 * L0)
            U[ln["burgers"]] += ln["sign"] * theta * bmag / (2.0 * np.pi)
            net[ln["burgers"]] += ln["sign"]
        if any(abs(v) > 1e-9 for v in net.values()):
            raise ValueError(f"net Burgers per component must cancel: {net}")
        self.psi = self._bcc_field(X - U["x"], Y - U["y"], Z - U["z"], amp=amp)
        self._fix_mean()

    def init_bicrystal(self, tilt_deg=15.0, amp=0.25):
        """Two BCC grains meeting at a flat (001) boundary at z=lz/2; the
        upper grain is rotated by `tilt_deg` about z. A clean single grain
        boundary lets the M12 detector resolve individual GB-emitted
        dislocation lines under tension (unlike a quenched polycrystal whose
        GB network percolates). Grains are large (= box halves), so lines
        are spatially isolated."""
        x = np.arange(self.nx) * self.dx
        y = np.arange(self.ny) * self.dy
        z = np.arange(self.nz) * self.dz
        Z, Y, X = np.meshgrid(z, y, x, indexing="ij")
        qx, qy, qz = self._snapped_q()
        th = np.deg2rad(tilt_deg)
        c, s = np.cos(th), np.sin(th)

        def bcc(Xr, Yr, Zr):
            return (np.cos(qx * Xr) * np.cos(qy * Yr)
                    + np.cos(qy * Yr) * np.cos(qz * Zr)
                    + np.cos(qz * Zr) * np.cos(qx * Xr))

        f_lo = bcc(X, Y, Z)
        Xr = c * X - s * Y
        Yr = s * X + c * Y
        f_hi = bcc(Xr, Yr, Z)
        xi = 1.5 * A_BCC
        w = 0.5 * (1.0 + np.tanh((Z - self.lz / 2.0) / xi))  # 0 lower,1 upper
        self.psi = self.psi_bar + amp * ((1 - w) * f_lo + w * f_hi)
        self._fix_mean()

    def init_bicrystal_csl(self, amp=0.25):
        """Σ5 twist bicrystal: lower grain unrotated, upper grain twisted about
        z by θ = atan(3/4) ≈ 36.87°. With exactly 5 BCC cells in x and y, the
        rotated wave vectors map onto the FFT grid (3,4,5 Pythagorean triple),
        so BOTH grains are box-commensurate and the only defects are at the
        single (001) twist boundary — no spurious bulk mismatch (cf. the
        non-CSL init_bicrystal). Requires nx == ny with 5 cells in-plane;
        caller must build the box with that commensurate dx (see
        run_bicrystal_csl)."""
        x = np.arange(self.nx) * self.dx
        y = np.arange(self.ny) * self.dy
        z = np.arange(self.nz) * self.dz
        Z, Y, X = np.meshgrid(z, y, x, indexing="ij")
        qx, qy, qz = self._snapped_q()
        c, s = 4.0 / 5.0, 3.0 / 5.0          # cos, sin of the Σ5 twist

        def bcc(Xr, Yr, Zr):
            return (np.cos(qx * Xr) * np.cos(qy * Yr)
                    + np.cos(qy * Yr) * np.cos(qz * Zr)
                    + np.cos(qz * Zr) * np.cos(qx * Xr))

        f_lo = bcc(X, Y, Z)
        f_hi = bcc(c * X - s * Y, s * X + c * Y, Z)
        xi = 1.5 * A_BCC
        w = 0.5 * (1.0 + np.tanh((Z - self.lz / 2.0) / xi))
        self.psi = self.psi_bar + amp * ((1 - w) * f_lo + w * f_hi)
        self._fix_mean()

    def _fix_mean(self):
        self.psi += self.psi_bar - self.psi.mean()

    # ---------- evolution ----------
    def step(self, dt, n=1):
        shape = self.psi.shape
        psi = self.psi
        # real field is the state; full round trip per iteration projects out
        # phantom spectral components (see pfc2d.step for the failure mode)
        for _ in range(n):
            psi_h = _rfftn(psi)
            # stabilized splitting (see pfc2d.step): C >= max 3 psi^2
            C = 3.0 * float(np.max(psi * psi))
            nl_h = _rfftn(psi ** 3)
            psi_h = ((psi_h - dt * self.k2 * (nl_h - C * psi_h))
                     / (1.0 + dt * self.k2 * (self.lin + C)))
            psi = _irfftn(psi_h, shape)
            self.time += dt
        self.psi = psi

    def step_mpfc(self, dt, n=1, beta=10.0):
        """Modified-PFC (Stefanovic inertial term) for fast elastic relaxation;
        beta=0 reduces exactly to step(). See pfc2d.step_mpfc."""
        shape = self.psi.shape
        psi = self.psi
        psi_prev = getattr(self, "_psi_prev", None)
        if psi_prev is None or psi_prev.shape != shape:
            psi_prev = psi.copy()
        a = beta / dt ** 2
        b = 1.0 / dt
        for _ in range(n):
            C = 3.0 * float(np.max(psi * psi))
            nl_h = _rfftn(psi ** 3)
            psi_h = _rfftn(psi)
            psi_prev_h = _rfftn(psi_prev)
            rhs = (a * (2.0 * psi_h - psi_prev_h) + b * psi_h
                   - self.k2 * (nl_h - C * psi_h))
            psi_new_h = rhs / (a + b + self.k2 * (self.lin + C))
            psi_prev = psi
            psi = _irfftn(psi_new_h, shape)
            self.time += dt
        self.psi = psi
        self._psi_prev = psi_prev

    def free_energy(self):
        psi_h = _rfftn(self.psi)
        lin_term = _irfftn(self.lin * psi_h, self.psi.shape)
        f = 0.5 * self.psi * lin_term + 0.25 * self.psi ** 4
        return float(f.mean())

    def stress(self, deps=1e-4):
        """d f / d exx along the volume-conserving tension path."""
        e0 = (self.exx, self.eyy, self.ezz)
        try:
            self.exx = (1.0 + e0[0]) * (1.0 + deps) - 1.0
            self.eyy = self.ezz = 1.0 / np.sqrt(1.0 + self.exx) - 1.0
            self._update_k()
            fp = self.free_energy()
            self.exx = (1.0 + e0[0]) * (1.0 - deps) - 1.0
            self.eyy = self.ezz = 1.0 / np.sqrt(1.0 + self.exx) - 1.0
            self._update_k()
            fm = self.free_energy()
        finally:
            self.exx, self.eyy, self.ezz = e0
            self._update_k()
        return (fp - fm) / (2.0 * deps * (1.0 + e0[0]))

    # ---------- io ----------
    def save(self, path):
        np.savez_compressed(
            path, psi=self.psi, r=self.r, psi_bar=self.psi_bar,
            dx0=self.dx0, exx=self.exx, eyy=self.eyy, ezz=self.ezz,
            gxz=self.gxz, time=self.time)

    @classmethod
    def load(cls, path):
        d = np.load(path)
        nz, ny, nx = d["psi"].shape
        m = cls(nx, ny, nz, dx=float(d["dx0"]), r=float(d["r"]),
                psi_bar=float(d["psi_bar"]))
        m.psi = d["psi"]
        m.exx, m.eyy, m.ezz = (float(d["exx"]), float(d["eyy"]),
                               float(d["ezz"]))
        m.gxz = float(d["gxz"]) if "gxz" in d else 0.0
        m.time = float(d["time"])
        m._update_k()
        return m


def find_peaks_3d(psi, dx, dy, dz, a0=A_BCC, refine=2):
    """Atom positions from 3D density maxima, with sub-grid refinement:
    windowed center-of-mass around each maximum (grid-quantized positions
    bias nearest-neighbor distances ~5-7% low)."""
    from scipy.ndimage import maximum_filter, label, center_of_mass
    size = max(3, int(0.4 * a0 / max(dx, dy, dz)))
    mx = maximum_filter(psi, size=size, mode="wrap")
    mask = (psi >= mx) & (psi > psi.mean() + 0.3 * psi.std())
    lbl, n = label(mask)
    if n == 0:
        return np.zeros((0, 3))
    coms = center_of_mass(psi - psi.min(), lbl, range(1, n + 1))
    nz, ny, nx = psi.shape
    w = refine
    pts = []
    base = psi.min()
    for c in coms:
        iz, iy, ix = (int(round(c[0])) % nz, int(round(c[1])) % ny,
                      int(round(c[2])) % nx)
        zi = np.arange(iz - w, iz + w + 1)
        yi = np.arange(iy - w, iy + w + 1)
        xi = np.arange(ix - w, ix + w + 1)
        block = psi[np.ix_(zi % nz, yi % ny, xi % nx)] - base
        wsum = block.sum()
        Z, Y, X = np.meshgrid(zi, yi, xi, indexing="ij")
        pts.append([(X * block).sum() / wsum * dx,
                    (Y * block).sum() / wsum * dy,
                    (Z * block).sum() / wsum * dz])
    p = np.array(pts)
    for ax, L in enumerate((nx * dx, ny * dy, nz * dz)):
        p[:, ax] %= L
        # float modulo of a tiny negative can return exactly L
        p[:, ax] = np.where(p[:, ax] >= L, p[:, ax] - L, p[:, ax])
    return p
