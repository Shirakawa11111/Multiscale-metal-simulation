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
        self.k2 = KX ** 2 + KY ** 2 + KZ ** 2
        self.lin = self.r + (1.0 - self.k2) ** 2

    def apply_strain(self, dexx, volume_conserving=True):
        new_fx = (1.0 + self.exx) * (1.0 + dexx)
        self.exx = new_fx - 1.0
        if volume_conserving:
            f_t = 1.0 / np.sqrt(1.0 + self.exx)
            self.eyy = self.ezz = f_t - 1.0
        self._update_k()

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
            time=self.time)

    @classmethod
    def load(cls, path):
        d = np.load(path)
        nz, ny, nx = d["psi"].shape
        m = cls(nx, ny, nz, dx=float(d["dx0"]), r=float(d["r"]),
                psi_bar=float(d["psi_bar"]))
        m.psi = d["psi"]
        m.exx, m.eyy, m.ezz = (float(d["exx"]), float(d["eyy"]),
                               float(d["ezz"]))
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
