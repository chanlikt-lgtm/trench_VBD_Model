#!/usr/bin/env python3
"""Consequence of the factor-2 stiffness error.

Re-solve the P2 linear problem with K/2 (the correct quadrature) and rerun the
Bernstein census. If the earlier P1/P2 activation discrepancy was caused by the
doubled operator, the corrected P2 solution should now show an overshoot.
"""
from __future__ import annotations
import json
import numpy as np
import scipy.sparse.linalg as spla

import p2_audit as A
from p2_bernstein import B

M = A.M
V_D = A.V_DRAIN


def run(fn, tag):
    z = np.load(fn)
    nodes, tris, mat = z["nodes"], z["tris"], z["mat"]
    x = nodes[:, 0]
    K, f = A.assemble(nodes, tris, mat)
    K = 0.5*K                      # correct the doubled quadrature weight

    si_all = mat == M.SI
    si_nodes = np.unique(tris[si_all])
    gate_nodes = np.unique(tris[mat == M.GATE])
    shield_nodes = np.unique(tris[mat == M.SHIELD])
    src = si_nodes[x[si_nodes] <= 1e-12]
    drn = si_nodes[x[si_nodes] >= M.BOTTOM_DEPTH - 1e-9]
    body = si_nodes[x[si_nodes] <= A.BODY_BOTTOM + 1e-12]

    fixed = np.zeros(len(nodes), bool); val = np.zeros(len(nodes))
    for s, v in ((src, A.V_SRC), (body, A.V_SRC), (drn, V_D),
                 (gate_nodes, A.V_G), (shield_nodes, A.V_SG)):
        fixed[s] = True; val[s] = v
    idx = np.flatnonzero(~fixed)
    phi = val.copy()
    phi[idx] = spla.spsolve(K[idx][:, idx].tocsc(), (f - K @ val)[idx])
    r = K @ phi - f; r[fixed] = 0.0
    R = float(np.linalg.norm(r[idx]) /
              (np.linalg.norm(K[idx][:, idx] @ phi[idx]) + np.linalg.norm(f[idx])))

    u = V_D - phi
    b = u[tris[si_all]] @ B.T
    free_si = si_nodes[~fixed[si_nodes]]
    imin = free_si[np.argmin(u[free_si])]
    print(f"\n[{tag}] corrected K/2")
    print(f"  R_free={R:.3e}  phi=[{phi.min():.6f},{phi.max():.6f}]")
    print(f"  min b = {b.min():+.6e}  negatives={int((b < -1e-12).sum())}  "
          f"elements affected={int((b < -1e-12).any(axis=1).sum())}/{int(si_all.sum())}")
    print(f"  min u free Si = {u[free_si].min():+.9f} V at "
          f"(x={nodes[imin,0]:.5f}, y={nodes[imin,1]:.5f})")
    print(f"  max(phi)-V_D = {(phi[si_nodes]-V_D).max():+.9f}")
    return dict(R_free=R, min_b=float(b.min()), neg=int((b < -1e-12).sum()),
                min_u=float(u[free_si].min()), u_loc=[float(nodes[imin, 0]), float(nodes[imin, 1])],
                overshoot=float((phi[si_nodes]-V_D).max()))


if __name__ == "__main__":
    out = {}
    for tag, fn in (("r3", "p2_r3_solution.npz"), ("r4", "p2_r4_solution.npz")):
        out[tag] = run(fn, tag)
    json.dump(out, open("p2_corrected_census.json", "w"), indent=2)
    print("\nwrote p2_corrected_census.json")
