"""Stage-2b: two-mode PFC (crystallography fix supplying dissociable
dislocations + sessile junction products). First milestone = the existence
gate: does the two-mode operator crystallize the target lattice?

2D square lattice testbed (cheapest two-mode case; intersecting slip systems +
richer junction geometry than one-mode triangular). Reciprocal families at
k=1 ({10}) and k=sqrt(2) ({11}). The linear operator
   L(k) = r + (k^2-1)^2 (k^2-2)^2
has degenerate wells at k=1 and k=sqrt(2) (curvature ~matches one-mode at k=1),
liquid-stable at k=0 (L=4>0). Climb-suppressed mobility (Stage-2a) composes via
glide_kc, so this module can test the DUAL fix (two-mode + climb suppression).
"""

import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
from pfc2d import PFC2D, _rfft2, _irfft2


class TwoModePFC(PFC2D):
    # robust square-lattice regime from the (r,psi_bar) stability scan
    # (gap=0.023, amp=1.31, melt-quench crystallizes): r=-0.45, psi_bar=-0.30
    def __init__(self, nx, ny, dx=np.pi / 4, r=-0.45, psi_bar=-0.30):
        super().__init__(nx, ny, dx=dx, r=r, psi_bar=psi_bar)

    def square_field(self, X, Y, amp=0.35, amp2=0.18, ux=None, uy=None):
        q = 2 * np.pi * max(1, round(1.0 * self.lx / (2 * np.pi))) / self.lx
        Xe = X - (ux if ux is not None else 0.0)
        Ye = Y - (uy if uy is not None else 0.0)
        f1 = np.cos(q * Xe) + np.cos(q * Ye)
        f2 = np.cos(q * (Xe + Ye)) + np.cos(q * (Xe - Ye))
        return self.psi_bar + amp * f1 + amp2 * f2

    def init_dislocations_square(self, cores, images=1):
        """Phase-winding dislocations in the two-mode SQUARE lattice. cores:
        (x_frac, y_frac, sign[, burgers_deg]); square slip systems are at
        0/90 deg (<10>) and 45/135 deg (<11>) — richer junction geometry than
        one-mode triangular. Net Burgers vector must cancel."""
        norm = [(c[0], c[1], c[2], c[3] if len(c) > 3 else 0.0) for c in cores]
        bx = sum(s * np.cos(np.radians(a)) for _, _, s, a in norm)
        by = sum(s * np.sin(np.radians(a)) for _, _, s, a in norm)
        if abs(bx) > 1e-9 or abs(by) > 1e-9:
            raise ValueError("net Burgers vector must be zero")
        x = np.arange(self.nx) * self.dx
        y = np.arange(self.ny) * self.dy
        X, Y = np.meshgrid(x, y)
        q = 2 * np.pi * max(1, round(1.0 * self.lx / (2 * np.pi))) / self.lx
        b = 2 * np.pi / q
        ux = np.zeros_like(X)
        uy = np.zeros_like(Y)
        for sx in range(-images, images + 1):
            for sy in range(-images, images + 1):
                for fx, fy, s, adeg in norm:
                    th = np.arctan2(Y - fy * self.ly + sy * self.ly,
                                    X - fx * self.lx + sx * self.lx)
                    a = np.radians(adeg)
                    ux += s * th * np.cos(a)
                    uy += s * th * np.sin(a)
        ux *= b / (2 * np.pi)
        uy *= b / (2 * np.pi)
        self.psi = self.square_field(X, Y, ux=ux, uy=uy)
        self._fix_mean()

    def _update_k(self):
        kx = 2.0 * np.pi * np.fft.rfftfreq(self.nx, d=self.dx)
        ky = 2.0 * np.pi * np.fft.fftfreq(self.ny, d=self.dy)
        KX, KY = np.meshgrid(kx, ky)
        kyp = KY - self.gamma * KX
        self.k2 = KX ** 2 + kyp ** 2
        # two-mode operator: wells at k^2=1 and k^2=2
        self.lin = self.r + (self.k2 - 1.0) ** 2 * (self.k2 - 2.0) ** 2
        if self.glide_kc > 0:
            self.mob = self.k2 ** 2 / (self.k2 + self.glide_kc ** 2)
        else:
            self.mob = self.k2

    def init_square(self, amp=0.2, amp2=0.1, noise=0.0, seed=0):
        """Seed a two-mode square lattice: first mode {10} at k=1, second mode
        {11} at k=sqrt(2). Wave numbers snapped to the box."""
        x = np.arange(self.nx) * self.dx
        y = np.arange(self.ny) * self.dy
        X, Y = np.meshgrid(x, y)
        q = 2 * np.pi * max(1, round(1.0 * self.lx / (2 * np.pi))) / self.lx
        f1 = np.cos(q * X) + np.cos(q * Y)
        f2 = np.cos(q * (X + Y)) + np.cos(q * (X - Y))
        self.psi = self.psi_bar + amp * f1 + amp2 * f2
        if noise > 0:
            self.psi += noise * np.random.default_rng(seed).standard_normal(
                self.psi.shape)
        self._fix_mean()


if __name__ == "__main__":
    import time
    from defect_analysis import find_peaks, lattice_spacing
    t0 = time.time()
    # existence gate: seeded square crystal relaxes & stays square; energy below
    # liquid; melt-quench forms a lattice.
    n = 128
    m = TwoModePFC(n, n, r=-0.20, psi_bar=-0.30)
    f_liq = m.free_energy()
    m.init_square()
    e0 = m.free_energy()
    m.step(0.5, n=600)
    e1 = m.free_energy()
    pts = find_peaks(m.psi, m.dx, m.dy)
    # square lattice: nearest-neighbour spacing = 2pi/q ; check peak count
    q = 2 * np.pi * max(1, round(m.lx / (2 * np.pi))) / m.lx
    a_sq = 2 * np.pi / q
    n_expect = (m.lx * m.ly) / a_sq ** 2          # 1 atom per a_sq^2 (square)
    print(f"two-mode square: F liquid={f_liq:.5f} seed={e0:.5f} -> relaxed={e1:.5f}")
    print(f"  peaks={len(pts)} expect~{n_expect:.0f} (a_sq={a_sq:.3f}), "
          f"crystalline_below_liquid={e1 < f_liq}")
    # melt quench
    mq = TwoModePFC(n, n, r=-0.20, psi_bar=-0.30)
    mq.init_random(noise=0.05, seed=2)
    mq.step(0.5, n=2000)
    eq = mq.free_energy()
    ptsq = find_peaks(mq.psi, mq.dx, mq.dy)
    print(f"  melt-quench: F={eq:.5f} peaks={len(ptsq)} "
          f"crystallized={eq < f_liq and len(ptsq) > 0.5 * n_expect}")
    ok = (e1 < f_liq and len(pts) > 0.7 * n_expect
          and abs(np.abs(m.psi - m.psi.mean()).max()) > 0.1)
    print(f"wall {time.time()-t0:.0f}s  EXISTENCE GATE: {'PASS' if ok else 'FAIL'}")
