#!/usr/bin/env python3
"""Bernstein census on both P2 resolutions, including the margin to activation
measured on FREE nodes only (Dirichlet nodes are at u=0 by construction)."""
from __future__ import annotations
import numpy as np
import p2_audit as A
from p2_bernstein import B

M = A.M
V_D = A.V_DRAIN


def census(fn, tag):
    z = np.load(fn)
    nodes, tris, mat, phi, fixed = z["nodes"], z["tris"], z["mat"], z["phi"], z["fixed"]
    u = V_D - phi
    si = np.flatnonzero(mat == M.SI)
    tsi = tris[si]
    b = u[tsi] @ B.T
    neg = b < -1e-12
    free = ~fixed
    si_nodes = np.unique(tsi)
    free_si = si_nodes[free[si_nodes]]
    print(f"[{tag}] Si elements={len(si)} free Si nodes={len(free_si)}")
    print(f"  min b = {b.min():+.6e}   negatives={int(neg.sum())}   "
          f"elements with any negative={int(neg.any(axis=1).sum())}")
    print(f"  min u over free Si nodes = {u[free_si].min():+.6e} V "
          f"(x={nodes[free_si[np.argmin(u[free_si])],0]:.4f}, "
          f"y={nodes[free_si[np.argmin(u[free_si])],1]:.4f})")
    print(f"  number of free Si nodes with u < 0.5 V: {int((u[free_si] < 0.5).sum())}")
    print(f"  number of free Si nodes with u < 1.0 V: {int((u[free_si] < 1.0).sum())}")
    return dict(min_b=float(b.min()), neg=int(neg.sum()),
                min_u_free=float(u[free_si].min()))


if __name__ == "__main__":
    out = {}
    for tag, fn in (("r3", "p2_r3_solution.npz"), ("r4", "p2_r4_solution.npz")):
        out[tag] = census(fn, tag)
    import json
    json.dump(out, open("p2_activation_margin.json", "w"), indent=2)
    print("wrote p2_activation_margin.json")
