#!/usr/bin/env python3
"""Manufactured active-obstacle test, solved with a primal active-set QP.

trust-constr failed to converge on this problem, so this replaces it with the
standard primal active-set method for a convex QP:

    min 1/2 x'Hx - g'x   s.t.  Ax >= b

Start from a Phase-I feasible point (HiGHS LP), maintain a working set W,
solve the equality-constrained KKT system for the step, and drop/add
constraints by multiplier sign and blocking ratio.
"""
from __future__ import annotations
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import linprog

import p2_audit as A
from p2_bernstein import B
from p2_obstacle_test import build_mesh, assemble

TOL = 1e-9


def active_set_qp(H, g, Aineq, bineq, x0, max_iter=400, verbose=True):
    n = H.shape[0]
    x = x0.copy()
    W = [i for i in range(len(bineq)) if (Aineq[i] @ x - bineq[i]) < 1e-8]
    for it in range(max_iter):
        AW = Aineq[W] if W else np.zeros((0, n))
        KKT = np.block([[H, -AW.T], [AW, np.zeros((len(W), len(W)))]])
        rhs = np.concatenate([g - H @ x, np.zeros(len(W))])
        sol, *_ = np.linalg.lstsq(KKT, rhs, rcond=None)
        p, lam = sol[:n], sol[n:]
        if np.linalg.norm(p) < 1e-11:
            if len(W) == 0 or lam.min() >= -1e-9:
                if verbose:
                    print(f"    converged in {it} iterations, |W|={len(W)}")
                return x, W, lam
            j = int(np.argmin(lam))
            W.pop(j)
            continue
        Ap = Aineq @ p
        slack = Aineq @ x - bineq
        block = Ap < -1e-14
        if block.any():
            ratios = slack[block] / (-Ap[block])
            k = int(np.argmin(ratios))
            alpha = min(1.0, float(ratios[k]))
        else:
            alpha, k = 1.0, None
        x = x + alpha * p
        if alpha < 1.0 - 1e-13 and k is not None:
            ci = int(np.flatnonzero(block)[k])
            if ci not in W:
                W.append(ci)
    if verbose:
        print(f"    WARNING: max iterations, |W|={len(W)}")
    return x, W, lam


def main():
    nodes, tris, mat = build_mesh()
    n = len(nodes)
    xcoord = nodes[:, 0]
    K = assemble(nodes, tris, mat)

    fixed = np.zeros(n, bool); val = np.zeros(n)
    left = xcoord <= 1e-12; right = xcoord >= 1-1e-12
    fixed[left] = True; val[left] = 1.0
    fixed[right] = True; val[right] = 0.0
    idx = np.flatnonzero(~fixed)

    # Negative source in silicon drives the interior below zero while every
    # Dirichlet value stays nonnegative (a negative Dirichlet value would make
    # the vertex Bernstein row b_v = u_v infeasible by construction).
    _, _, area, _, _, _ = A.element_geom(nodes, tris)
    f = np.zeros(n)
    S = 200.0
    for e in np.flatnonzero(mat == A.M.SI):
        f[tris[e][3:]] += -S * abs(area[e]) / 3.0

    u = val.copy()
    u[idx] = spla.spsolve(K[idx][:, idx].tocsc(),
                          f[idx] - K[idx][:, fixed] @ val[fixed])
    si = np.flatnonzero(mat == A.M.SI)
    bu = u[tris[si]] @ B.T
    print(f"[U] unconstrained: min b = {bu.min():+.6e}  negatives = {int((bu < 0).sum())} "
          f"elements affected = {int((bu < 0).any(axis=1).sum())}/{len(si)}")

    C = np.zeros((6*len(si), n))
    for e, el in enumerate(si):
        C[6*e:6*e+6, tris[el]] = B
    C = sp.csr_matrix(C)
    C_F = C[:, idx]
    cD = C[:, fixed] @ val[fixed]

    H = K[idx][:, idx].toarray()
    g = f[idx] - K[idx][:, fixed] @ val[fixed]
    Aineq = C_F.toarray()
    bineq = -cD

    lp = linprog(np.zeros(len(idx)), A_ub=-Aineq, b_ub=cD,
                 bounds=(None, None), method="highs")
    print(f"[phase-I] {lp.status}: {lp.message}")
    x0 = lp.x

    print("[QP] active-set")
    z, W, lam = active_set_qp(H, g, Aineq, bineq, x0)
    uq = val.copy(); uq[idx] = z
    slack = Aineq @ z - bineq
    print(f"[P] min slack = {slack.min():+.3e}   active = {int((slack < 1e-7).sum())}"
          f" / {len(slack)}")
    print(f"    u range [{uq.min():+.6e}, {uq.max():+.6f}]")

    r = H @ z - g
    lam_full = np.zeros(len(slack)); lam_full[W] = lam
    stat = r - Aineq.T @ lam_full
    print(f"[S] ||stationarity|| = {np.linalg.norm(stat):.3e}   "
          f"relative = {np.linalg.norm(stat)/max(np.linalg.norm(r),1e-30):.3e}")
    print(f"[D] lam_min = {lam.min():+.3e}   lam_max = {lam.max():.3e}")
    comp = np.max(np.abs(lam_full * slack))
    print(f"[C] complementarity max|lam*slack| = {comp:.3e}")

    ox_nodes = np.unique(tris[mat == A.M.OX])
    ox_free = np.intersect1d(ox_nodes, idx)
    r_all = K @ uq
    print(f"[O] oxide free rows max|K u| = {np.abs(r_all[ox_free]).max():.3e}")

    worst = 1e9
    L1s, L2s = np.meshgrid(np.linspace(0.01, 0.97, 14), np.linspace(0.01, 0.97, 14))
    keep = (L1s + L2s) < 0.99
    P1, P2 = L1s[keep], L2s[keep]
    for el in si:
        vals = uq[tris[el]]
        for a, b in zip(P1, P2):
            worst = min(worst, float(np.dot(vals, A.shape_vals(a, b))))
    print(f"[I] dense sampling of {len(si)} Si elements: min u = {worst:+.3e}")

    np.savez("p2_obstacle_test.npz", nodes=nodes, tris=tris, mat=mat,
             u_uncon=u, u_con=uq, slack=slack, W=np.array(W))


if __name__ == "__main__":
    main()
