#!/usr/bin/env python3
"""Manufactured active-obstacle test solved with OSQP.

min 1/2 x'Hx - g'x   s.t.  C_F x + c_D >= 0

Same two-material patch as p2_obstacle_test: oxide left (unconstrained),
silicon right (Bernstein b >= 0), nonnegative Dirichlet data, negative silicon
source so the unconstrained solution violates the constraint naturally.
"""
from __future__ import annotations
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import osqp

import p2_audit as A
from p2_bernstein import B
from p2_obstacle_test import build_mesh, assemble


def main():
    nodes, tris, mat = build_mesh()
    n = len(nodes)
    xc = nodes[:, 0]
    K = assemble(nodes, tris, mat)

    fixed = np.zeros(n, bool); val = np.zeros(n)
    fixed[xc <= 1e-12] = True; val[xc <= 1e-12] = 1.0
    fixed[xc >= 1-1e-12] = True; val[xc >= 1-1e-12] = 0.0
    idx = np.flatnonzero(~fixed)

    _, _, area, _, _, _ = A.element_geom(nodes, tris)
    f = np.zeros(n)
    S = 60.0
    for e in np.flatnonzero(mat == A.M.SI):
        f[tris[e][3:]] += -S * abs(area[e]) / 3.0

    u = val.copy()
    u[idx] = spla.spsolve(K[idx][:, idx].tocsc(), f[idx] - K[idx][:, fixed] @ val[fixed])

    si = np.flatnonzero(mat == A.M.SI)
    bu = u[tris[si]] @ B.T
    print(f"[U] unconstrained: min b = {bu.min():+.6e}  negatives = {int((bu<0).sum())} "
          f"elements affected = {int((bu<0).any(axis=1).sum())}/{len(si)}")
    print(f"    u range [{u.min():+.6f}, {u.max():+.6f}]")

    C = np.zeros((6*len(si), n))
    for e, el in enumerate(si):
        C[6*e:6*e+6, tris[el]] = B
    C = sp.csr_matrix(C)
    C_F = C[:, idx].tocsc()
    cD = C[:, fixed] @ val[fixed]

    H = sp.csc_matrix(K[idx][:, idx])
    g = f[idx] - K[idx][:, fixed] @ val[fixed]

    prob = osqp.OSQP()
    prob.setup(P=H, q=-g, A=C_F, l=-cD, u=np.full(len(cD), np.inf),
               verbose=False, eps_abs=1e-11, eps_rel=1e-11,
               max_iter=200000, polishing=True)
    res = prob.solve()
    print(f"\n[QP] {res.info.status}  iters={res.info.iter}  "
          f"obj={res.info.obj_val:.10e}")
    z = res.x
    y = res.y

    uq = val.copy(); uq[idx] = z
    slack = C_F @ z + cD
    print(f"[P] min slack = {slack.min():+.3e}   active (<1e-7) = {int((slack<1e-7).sum())}"
          f" / {len(slack)}")
    print(f"    u range [{uq.min():+.6e}, {uq.max():+.6f}]")

    r = H @ z - g
    for name, stat in (("P x + q + A'y", r + C_F.T @ y),
                       ("P x + q - A'y", r - C_F.T @ y)):
        print(f"[S] ||{name}|| / ||r|| = {np.linalg.norm(stat)/np.linalg.norm(r):.3e}")
    lam_osqp = -y
    print(f"[D] lambda = -y: min = {lam_osqp.min():+.3e}  max = {lam_osqp.max():.3e}  "
          f"negative count = {int((lam_osqp < -1e-9).sum())}")
    comp = np.max(np.abs(y * slack))
    print(f"[C] complementarity max|y*slack| = {comp:.3e}")

    si_nodes = set(np.unique(tris[si]).tolist())
    ox_nodes = [i for i in np.unique(tris[mat == A.M.OX]) if int(i) not in si_nodes]
    ox_free = np.intersect1d(np.array(ox_nodes, dtype=int), idx)
    r_all = K @ uq
    print(f"[O] pure-oxide free rows max|K u| = {np.abs(r_all[ox_free]).max():.3e} "
          f"({len(ox_free)} nodes)")

    worst = 1e9
    L1s, L2s = np.meshgrid(np.linspace(0.005, 0.995, 16), np.linspace(0.005, 0.995, 16))
    keep = (L1s + L2s) < 0.995
    P1, P2 = L1s[keep], L2s[keep]
    for el in si:
        vals = uq[tris[el]]
        for a, b in zip(P1, P2):
            worst = min(worst, float(np.dot(vals, A.shape_vals(a, b))))
    print(f"[I] dense sampling of {len(si)} Si elements: min u = {worst:+.3e}")

    out = dict(status=str(res.info.status), iters=int(res.info.iter),
               unconstrained_min_b=float(bu.min()),
               min_slack=float(slack.min()),
               active=int((slack < 1e-7).sum()),
               lam_min=float(lam_osqp.min()), comp=float(comp),
               oxide_res=float(np.abs(r_all[ox_free]).max()),
               sample_min=float(worst))
    import json
    json.dump(out, open("p2_obstacle_test.json", "w"), indent=2)
    print("wrote p2_obstacle_test.json")


if __name__ == "__main__":
    main()
