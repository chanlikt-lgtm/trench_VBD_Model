#!/usr/bin/env python3
"""F2: provable lower/upper bounds for the P2 neutral area.

For each charged triangle, form the Bernstein coefficients of (u - eps).
Convex-hull property:
  all <= 0  -> entirely neutral  (add area to LOWER)
  all >  0  -> entirely outside  (add nothing)
  mixed     -> subdivide (de Casteljau at midpoints) until the unresolved
               area is below tol; unresolved area is added to UPPER only.
Then LOWER <= A_true <= UPPER.
"""
from __future__ import annotations
import json
import numpy as np

import p2_audit as A
from p2_bernstein import B

M = A.M
EPS = 1e-8
AREA_TOL = 1e-11
MAX_DEPTH = 24


def eval_bern(c, L):
    L1, L2, L3 = L
    Bs = (L1*L1, L2*L2, L3*L3, 2*L1*L2, 2*L2*L3, 2*L3*L1)
    return float(sum(ci*bi for ci, bi in zip(c, Bs)))


def sub_coeff(c, A_, B_, C_):
    def mid(P, Q):
        return ((P[0]+Q[0])/2, (P[1]+Q[1])/2, (P[2]+Q[2])/2)
    M01 = mid(A_, B_); M12 = mid(B_, C_); M20 = mid(C_, A_)
    pts = [A_, B_, C_, M01, M12, M20]
    v = [eval_bern(c, p) for p in pts]
    return [v[0], v[1], v[2],
            2*v[3] - (v[0]+v[1])/2,
            2*v[4] - (v[1]+v[2])/2,
            2*v[5] - (v[2]+v[0])/2]


def recurse(pc, A_, B_, C_, area, depth, acc):
    # coefficients of the CURRENT sub-triangle, from the parent polynomial
    # evaluated at this sub-triangle's nodes (all in parent barycentrics)
    sc = sub_coeff(pc, A_, B_, C_)
    if all(ci <= 0.0 for ci in sc):
        acc[0] += area; return
    if all(ci > 0.0 for ci in sc):
        return
    if area < AREA_TOL or depth >= MAX_DEPTH:
        acc[1] += area; return
    def mid(P, Q):
        return ((P[0]+Q[0])/2, (P[1]+Q[1])/2, (P[2]+Q[2])/2)
    M01 = mid(A_, B_); M12 = mid(B_, C_); M20 = mid(C_, A_)
    for P, Q, R in ((A_, M01, M20), (M01, B_, M12),
                    (M20, M12, C_), (M01, M12, M20)):
        recurse(pc, P, Q, R, area/4.0, depth+1, acc)


def bounds(nodes, tris, mat, u):
    area = np.abs(A.element_geom(nodes, tris)[2])
    cxs = nodes[tris][:, :3, 0].mean(axis=1)
    charged = (mat == M.SI) & (cxs > A.BODY_BOTTOM + 1e-12)
    els = np.flatnonzero(charged)
    lower = 0.0; unresolved = 0.0
    n_mixed = 0
    for e in els:
        c = B @ u[tris[e]] - EPS
        if not (all(ci <= 0.0 for ci in c) or all(ci > 0.0 for ci in c)):
            n_mixed += 1
        acc = [0.0, 0.0]
        recurse(c, (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0),
                area[e], 0, acc)
        lower += acc[0]
        unresolved += acc[1]
    upper = lower + unresolved
    return lower, upper, int(charged.sum()), n_mixed


def main():
    out = {}
    for tag, fn in (("r3", "p2_p3_r3_polished.npz"),
                    ("r4", "p2_p3_r4_polished.npz"),
                    ("mid", "p2_p3_rmid_polished.npz")):
        z = np.load(fn)
        nodes, tris, mat, u = z["nodes"], z["tris"], z["mat"], z["u"]
        lo, up, nch, nm = bounds(nodes, tris, mat, u)
        print(f"[{tag}] charged={nch} mixed={nm}  lower={lo:.12f}  upper={up:.12f}  "
              f"gap={up-lo:.3e}")
        out[tag] = dict(lower=lo, upper=up, gap=up-lo, charged=nch, mixed=nm)
    json.dump(out, open("p2_area_bounds.json", "w"), indent=2)
    print("\nwrote p2_area_bounds.json")


if __name__ == "__main__":
    main()
