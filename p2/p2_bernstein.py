#!/usr/bin/env python3
"""P2.2 — Bernstein positivity constraint layer (separately testable).

Transform (per silicon P2 triangle, node order v0,v1,v2,m01,m12,m20):
    b_v  = u_v
    b_m  = 2 u_m - (u_a + u_b)/2

Tests, in order:
  1. transform exactness on random quadratic polynomials
  2. sufficiency: b >= 0 implies u >= 0 inside (dense independent sampling)
  3. nodal non-negativity is NOT sufficient (explicit counterexample)
  4. Bernstein is sufficient but NOT necessary (explicit counterexample)
  5. constraint matrix on the real r3 silicon mesh: size, nnz
  6. violation census of the unconstrained linear P2 solution

u = V_D - phi, so u >= 0 enforces phi <= V_D.
"""
from __future__ import annotations
import json
import numpy as np
import scipy.sparse as sp

import p2_audit as A

M = A.M
V_D = A.V_DRAIN

# b = B u
B = np.array([
    [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
    [-0.5, -0.5, 0.0, 2.0, 0.0, 0.0],
    [0.0, -0.5, -0.5, 0.0, 2.0, 0.0],
    [-0.5, 0.0, -0.5, 0.0, 0.0, 2.0],
])


def bernstein_eval(b, L1, L2):
    """Evaluate the degree-2 Bernstein expansion at barycentric (L1,L2,L3)."""
    L3 = 1.0 - L1 - L2
    Bs = np.array([L1*L1, L2*L2, L3*L3, 2*L1*L2, 2*L2*L3, 2*L3*L1])
    return float(np.dot(b, Bs))


def nodal_from_poly(coef, XY):
    """Nodal values of a quadratic polynomial at six P2 nodes.

    coef = (c0, cx, cy, cxx, cxy, cyy);  p = c0 + cx x + cy y + cxx x^2 + cxy x y + cyy y^2
    XY: (6,2) node coordinates in the standard order (3 vertices then 3 midsides)
    """
    c0, cx, cy, cxx, cxy, cyy = coef
    x, y = XY[:, 0], XY[:, 1]
    return c0 + cx*x + cy*y + cxx*x*x + cxy*x*y + cyy*y*y


def sample_bary(n):
    """Random interior barycentric points (uniform on the triangle)."""
    r1 = np.sqrt(np.random.rand(n))
    r2 = np.random.rand(n)
    return r1*(1-r2), r1*r2


def main():
    rng = np.random.default_rng(20261019)
    out = {}

    # ---------- 1. transform exactness ----------
    errs = []
    for _ in range(200):
        P = rng.normal(size=(3, 2))
        XY = np.vstack([P, 0.5*(P[0]+P[1]), 0.5*(P[1]+P[2]), 0.5*(P[2]+P[0])])
        coef = rng.normal(size=6)
        u = nodal_from_poly(coef, XY)
        b = B @ u
        L1, L2 = sample_bary(200)
        ex = nodal_from_poly(coef, np.column_stack([
            L1*P[0, 0] + L2*P[1, 0] + (1-L1-L2)*P[2, 0],
            L1*P[0, 1] + L2*P[1, 1] + (1-L1-L2)*P[2, 1]]))
        ap = np.array([bernstein_eval(b, a, c) for a, c in zip(L1, L2)])
        errs.append(np.abs(ap - ex).max())
    out["transform_exactness"] = float(max(errs))
    print(f"[1] transform exactness: max |bernstein - exact quadratic| = {max(errs):.3e}")

    # ---------- 2. sufficiency ----------
    mins = []
    for _ in range(200):
        b = rng.random(6)
        L1, L2 = sample_bary(300)
        vals = [bernstein_eval(b, a, c) for a, c in zip(L1, L2)]
        mins.append(min(vals))
    out["sufficiency_min"] = float(min(mins))
    print(f"[2] sufficiency: min interior value over b>=0 samples = {min(mins):.3e} "
          f"(must be >= 0)")

    # ---------- 3. nodal non-negativity insufficient ----------
    u_bad = np.array([1.0, 1.0, 0.0, 0.0, 0.0, 0.0])
    b_bad = B @ u_bad
    L1, L2 = sample_bary(3000)
    vals = np.array([bernstein_eval(b_bad, a, c) for a, c in zip(L1, L2)])
    # u = N0 + N1 = L1(2L1-1) + L2(2L2-1); stationary at L1=L2=1/4 -> -1/4.
    print(f"[3] counterexample: u nodal all >=0, b = {np.array2string(b_bad, precision=2)}, "
          f"interior min = {vals.min():.6f}  (analytic min -1/4 = -0.25; "
          f"centroid -2/9 = {-2/9:.6f})")
    out["counter_nodal"] = dict(b=[float(v) for v in b_bad], min=float(vals.min()))

    # ---------- 4. Bernstein not necessary ----------
    b_nec = np.array([1.0, 1.0, 1.0, -1.0, 1.0, -1.0])
    vals = np.array([bernstein_eval(b_nec, a, c) for a, c in zip(L1, L2)])
    print(f"[4] counterexample: b has negatives {b_nec.tolist()}, interior min = "
          f"{vals.min():.3e} (function is (2L1-1)^2 >= 0)")
    out["counter_necessity"] = dict(b=[float(v) for v in b_nec], min=float(vals.min()))

    # ---------- 5. constraint matrix on the real mesh ----------
    z = np.load("p2_r3_solution.npz")
    nodes, tris, mat, phi = z["nodes"], z["tris"], z["mat"], z["phi"]
    si = np.flatnonzero(mat == M.SI)
    tsi = tris[si]
    n = len(nodes)
    rows = np.repeat(np.arange(6*len(si)), 6)
    cols = np.tile(tsi, (1, 6)).ravel()
    vals = np.tile(B.ravel(), len(si))
    C = sp.coo_matrix((vals, (rows, cols)), shape=(6*len(si), n)).tocsr()
    out["constraint_matrix"] = dict(rows=int(C.shape[0]), cols=int(C.shape[1]),
                                    nnz=int(C.nnz))
    print(f"[5] constraint matrix: {C.shape[0]} rows x {C.shape[1]} cols, nnz={C.nnz}")

    # ---------- 6. violation census of the unconstrained solution ----------
    u = V_D - phi
    ub = u[tsi]
    bs = ub @ B.T
    tol = 1e-9
    neg = bs < -tol
    n_elem_bad = int(neg.any(axis=1).sum())
    print(f"[6] u = V_D - phi : min u = {u.min():.6e}  max u = {u.max():.6f}")
    print(f"[6] Bernstein coefficients: min = {bs.min():.6e}  "
          f"negative entries = {int(neg.sum())}  elements with any negative = {n_elem_bad}"
          f" / {len(si)}")
    if neg.any():
        k = np.unravel_index(np.argmin(bs), bs.shape)
        el = si[k[0]]
        print(f"    worst: element {el} local coeff {k[1]} value {bs[k]:.6e} "
              f"centroid ({nodes[tsi[k[0]]][:3,0].mean():.4f}, {nodes[tsi[k[0]]][:3,1].mean():.4f})")
        # how negative, distribution
        vals_neg = bs[neg]
        print(f"    negative coefficient stats: min={vals_neg.min():.6e} "
              f"median={np.median(vals_neg):.6e} p95={np.percentile(vals_neg,5):.6e}")
    out["census"] = dict(min_u=float(u.min()), min_b=float(bs.min()),
                         n_negative=int(neg.sum()), n_elements=int(n_elem_bad),
                         n_si=int(len(si)))
    with open("p2_bernstein.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print("wrote p2_bernstein.json")


if __name__ == "__main__":
    main()
