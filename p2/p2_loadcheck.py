#!/usr/bin/env python3
"""Is the P2 charge term active, and which way does it push?

Solve the same system with and without the silicon charge load and compare.
If phi_load == phi_harmonic, the load is negligible. If the difference is
negative near the drain, the sign is opposite to the P1 convention.
"""
from __future__ import annotations
import numpy as np
import scipy.sparse.linalg as spla
import p2_audit as A

M = A.M
V_D = A.V_DRAIN


def main():
    z = np.load("p2_r3_solution.npz")
    nodes, tris, mat, phi = z["nodes"], z["tris"], z["mat"], z["phi"]
    x = nodes[:, 0]
    K, f = A.assemble(nodes, tris, mat)

    print(f"[load] sum(f) = {f.sum():.6e}   max(f) = {f.max():.6e}   "
          f"nnz(f) = {int((f != 0).sum())}")
    si = mat == M.SI
    area_si = float(np.sum(np.abs(A.element_geom(nodes, tris)[2])[si]))
    print(f"[load] Si area = {area_si:.6f} um^2   RHO = {A.RHO:.6f} V/um^2   "
          f"RHO*area = {A.RHO*area_si:.6f}  (sum f should match)")

    # Dirichlet set, as in p2_solve.py
    si_n = set(np.unique(tris[si]).tolist())
    ox_n = set(np.unique(tris[mat == M.OX]).tolist())
    gate_n = set(np.unique(tris[mat == M.GATE]).tolist())
    shield_n = set(np.unique(tris[mat == M.SHIELD]).tolist())
    src = {i for i in si_n if x[i] <= 1e-12}
    drn = {i for i in si_n if x[i] >= M.BOTTOM_DEPTH - 1e-9}
    body = {i for i in si_n if x[i] <= A.BODY_BOTTOM + 1e-12}
    fixed = np.zeros(len(nodes), bool); val = np.zeros(len(nodes))
    for s, v in ((src, A.V_SRC), (body, A.V_SRC), (drn, V_D),
                 (gate_n, A.V_G), (shield_n, A.V_SG)):
        for i in s:
            fixed[i] = True; val[i] = v
    idx = np.flatnonzero(~fixed)
    Kff = K[idx][:, idx].tocsc()

    def solve(rhs_full):
        p = val.copy()
        p[idx] = spla.spsolve(Kff, (rhs_full - K @ val)[idx])
        return p

    phi_load = phi                      # saved solution
    phi_zero = solve(np.zeros(len(nodes)))
    d = phi_load - phi_zero
    print(f"\n[compare] max|phi_load - phi_harmonic| = {np.abs(d).max():.6e} V")
    print(f"[compare] phi_loaded max = {phi_load.max():.6f}   "
          f"phi_harmonic max = {phi_zero.max():.6f}   V_D = {V_D:.6f}")
    for lo, hi, name in ((2.30, 2.42, "drain 2.30-2.42"), (1.5, 2.0, "mid 1.5-2.0"),
                         (0.6, 1.0, "upper 0.6-1.0")):
        m = (x > lo) & (x < hi)
        if m.any():
            print(f"  {name:16s} mean dphi = {d[m].mean():+.6e}  "
                  f"min u = {(V_D - phi_load[m]).min():+.6e}")

    u = V_D - phi_load
    print(f"\n[u] min={u.min():.6e}  argmin x={x[np.argmin(u)]:.4f}")
    print(f"[u] V_D - phi at drain nodes: "
          f"{np.array2string((V_D - phi_load[list(drn)])[:3], precision=6)}")

    # Is the interior maximum strictly below V_D, i.e. no overshoot at all?
    print(f"[overshoot] max(phi) - V_D = {phi_load.max() - V_D:+.6e}")


if __name__ == "__main__":
    main()
