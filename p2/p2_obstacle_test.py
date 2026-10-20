#!/usr/bin/env python3
"""P2.2 manufactured quadratic obstacle test with a NATURALLY active set.

Domain: unit square, P2 triangles, two materials split at x = 0.5
  oxide  x < 0.5   eps_r = 3.9   (unconstrained)
  Si     x > 0.5   eps_r = 11.7  (Bernstein b >= 0)

Dirichlet: u(0,y) = 1.0, u(1,y) = -0.2, top/bottom Neumann.
Interface value: 3.9(u-1) = 11.7(-0.2-u) -> u(0.5) = +0.1, so Si crosses zero
at x ~ 0.667 and part of the silicon is genuinely negative. The obstacle must
activate on its own; no node is forced.

Checks:
  U  unconstrained solution violates some Bernstein coefficients
  P  QP primal feasibility
  S  stationarity with nonnegative multipliers (recovered by NNLS on the active set)
  D  dual feasibility and complementarity
  O  oxide free rows satisfy unconstrained equilibrium
  I  dense interior sampling of every constrained silicon element is nonnegative
"""
from __future__ import annotations
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import minimize, LinearConstraint, nnls

import p2_audit as A
from p2_bernstein import B

EPS_SI, EPS_OX = A.EPS_SI, A.EPS_OX
QP, QW = A.QP, A.QW


def build_mesh(nx=8, ny=8, split=0.5):
    verts, vid = [], {}
    for j in range(ny+1):
        for i in range(nx+1):
            vid[(i, j)] = len(verts)
            verts.append((i/nx, j/ny))
    tris = []
    for j in range(ny):
        for i in range(nx):
            a, b = vid[(i, j)], vid[(i+1, j)]
            c, d = vid[(i+1, j+1)], vid[(i, j+1)]
            tris.append((a, b, c))
            tris.append((a, c, d))
    verts = np.array(verts, float)
    tris = np.array(tris)

    mid, mid_id = [], {}
    for t in tris:
        local = []
        for (p, q) in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            k = (min(p, q), max(p, q))
            if k not in mid_id:
                mid_id[k] = len(verts) + len(mid)
                mid.append(0.5*(verts[p] + verts[q]))
            local.append(mid_id[k])
        mid.append(None)  # placeholder guard, removed below
        mid.pop()
    nodes = np.vstack([verts, np.array(mid)])
    tris6 = np.empty((len(tris), 6), dtype=int)
    for e, t in enumerate(tris):
        ids = []
        for (p, q) in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            ids.append(mid_id[(min(p, q), max(p, q))])
        tris6[e] = [t[0], t[1], t[2], ids[0], ids[1], ids[2]]
    cx = nodes[tris6][:, :3, 0].mean(axis=1)
    mat = np.where(cx > split, A.M.SI, A.M.OX)
    return nodes, tris6, mat


def assemble(nodes, tris, mat):
    n = len(nodes)
    _, det, area, g1, g2, g3 = A.element_geom(nodes, tris)
    eps = np.where(mat == A.M.SI, EPS_SI, EPS_OX)
    rows = np.zeros((len(tris), 36), dtype=np.int64)
    cols = np.zeros((len(tris), 36), dtype=np.int64)
    vals = np.zeros((len(tris), 36))
    for m in range(len(tris)):
        ke = np.zeros((6, 6))
        for q in range(3):
            L1, L2 = QP[q]
            G = A.shape_grads(L1, L2, g1[m], g2[m], g3[m])
            ke += QW[q] * abs(det[m]) * eps[m] * (G @ G.T)
        ids = tris[m]
        rows[m] = np.repeat(ids, 6)
        cols[m] = np.tile(ids, 6)
        vals[m] = ke.ravel()
    K = sp.coo_matrix((vals.ravel(), (rows.ravel(), cols.ravel())),
                      shape=(n, n)).tocsr()
    return K


def main():
    nodes, tris, mat = build_mesh()
    n = len(nodes)
    x = nodes[:, 0]
    K = assemble(nodes, tris, mat)

    fixed = np.zeros(n, bool)
    val = np.zeros(n)
    left = x <= 1e-12
    right = x >= 1-1e-12
    fixed[left] = True; val[left] = 1.0
    fixed[right] = True; val[right] = -0.2
    free = ~fixed
    idx = np.flatnonzero(free)
    print(f"[mesh] nodes={n} elements={len(tris)} Si={int((mat==A.M.SI).sum())} "
          f"free={len(idx)} fixed={int(fixed.sum())}")

    # unconstrained
    u = val.copy()
    u[idx] = spla.spsolve(K[idx][:, idx].tocsc(), -(K[idx][:, fixed] @ val[fixed]))
    si = np.flatnonzero(mat == A.M.SI)
    bu = u[tris[si]] @ B.T
    print(f"[U] unconstrained: min b = {bu.min():+.6e}  negatives = {int((bu<0).sum())} "
          f"elements with any negative = {int((bu<0).any(axis=1).sum())}/{len(si)}")
    print(f"    unconstrained u range [{u.min():+.6f}, {u.max():+.6f}]")

    # constraints
    C = np.zeros((6*len(si), n))
    for e, el in enumerate(si):
        C[6*e:6*e+6, tris[el]] = B
    C = sp.csr_matrix(C)
    C_F = C[:, idx]; C_D = C[:, fixed]
    cD = C_D @ val[fixed]

    Kff = K[idx][:, idx].tocsc()
    Kfd = K[idx][:, fixed]
    g = -(Kfd @ val[fixed])          # objective 0.5 x'Kff x - g'x

    def fun(z):
        return 0.5*z @ (Kff @ z) - g @ z

    def jac(z):
        return Kff @ z - g

    def hess(z):
        return Kff

    res = minimize(fun, u[idx], jac=jac, hess=hess,
                   constraints=[LinearConstraint(C_F, 0.0, np.inf)],
                   method="trust-constr",
                   options=dict(gtol=1e-10, xtol=1e-12, maxiter=5000, verbose=0))
    z = res.x
    uq = val.copy(); uq[idx] = z
    slack = C_F @ z + cD
    print(f"\n[P] QP: {res.message}")
    print(f"    min slack = {slack.min():+.3e}   active (<1e-7) = {int((slack<1e-7).sum())}"
          f" / {len(slack)}")
    print(f"    constrained u range [{uq.min():+.6e}, {uq.max():+.6f}]")

    # stationarity with recovered multipliers on the active set
    r = Kff @ z - g
    act = slack < 1e-7
    lam, resid = nnls(C_F[act].T.toarray(), r)
    print(f"[S] stationarity: ||r||={np.linalg.norm(r):.3e}  "
          f"||C_act'lam - r||/||r||={resid/np.linalg.norm(r):.3e}  "
          f"lam_min={lam.min():+.3e}  lam_max={lam.max():.3e}")

    # complementarity
    comp = np.max(np.abs(lam * slack[act])) if act.any() else 0.0
    print(f"[D] dual min = {lam.min():+.3e}  complementarity max|lam*slack| = {comp:.3e}")

    # oxide equilibrium: free oxide rows carry no multiplier
    ox = np.flatnonzero(mat == A.M.SI)
    ox_free = np.flatnonzero(free & np.isin(np.arange(n), np.unique(tris[mat == A.M.OX])))
    r_all = K @ uq
    print(f"[O] oxide free rows: max|xK u| = {np.abs(r_all[ox_free]).max():.3e}  "
          f"(should be ~0)")

    # dense interior sampling
    worst = 1e9
    L1s, L2s = np.meshgrid(np.linspace(0.02, 0.96, 12), np.linspace(0.02, 0.96, 12))
    keep = (L1s + L2s) < 0.98
    P1, P2 = L1s[keep], L2s[keep]
    for el in si:
        ids = tris[el]
        vals = uq[ids]
        for a, b in zip(P1, P2):
            worst = min(worst, float(np.dot(vals, A.shape_vals(a, b))))
    print(f"[I] dense interior sampling of {len(si)} Si elements: min u = {worst:+.3e}")


if __name__ == "__main__":
    main()
