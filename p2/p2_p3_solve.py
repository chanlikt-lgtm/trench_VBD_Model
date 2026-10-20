#!/usr/bin/env python3
"""P2.3 — first physical trench obstacle solve.

Model (per the recorded decision):
  - corrected P2 operator (QW sums to 0.5)
  - P1-matched body: silicon elements with centroid x <= 0.54 carry no donor
    charge and every one of their nodes is pinned at u_D = V_D (equipotential)
  - gate and shield pinned at u_D = V_D; source pinned at V_D; whole drain
    plane pinned at u_D = 0
  - Bernstein b >= 0 on charged silicon elements only
  - variable u = V_D - phi, so u >= 0 enforces phi <= V_D

QP: min 1/2 u'K u - b'u  s.t. C u + c_D >= 0,  b = -f - K u_D
"""
from __future__ import annotations
import json
import numpy as np
import scipy.sparse as sp
import osqp

import p2_audit as A
from p2_bernstein import B

M = A.M
V_D = A.V_DRAIN


def main():
    z = np.load("p2_r3_fixed.npz")
    nodes, tris, mat = z["nodes"], z["tris"], z["mat"]
    x = nodes[:, 0]
    K, _ = A.assemble(nodes, tris, mat)
    _, _, area, _, _, _ = A.element_geom(nodes, tris)
    area = np.abs(area)          # element_geom returns signed area

    cxs = nodes[tris][:, :3, 0].mean(axis=1)
    si = mat == M.SI
    body = si & (cxs <= A.BODY_BOTTOM + 1e-12)
    charged = si & ~body

    f = np.zeros(len(nodes))
    for e in np.flatnonzero(charged):
        f[tris[e][3:]] += A.RHO*area[e]/3.0

    fixed = np.zeros(len(nodes), bool)
    val = np.zeros(len(nodes))
    body_nodes = np.unique(tris[body])
    si_nodes = np.unique(tris[si])
    gate = np.unique(tris[mat == M.GATE])
    shield = np.unique(tris[mat == M.SHIELD])
    src = si_nodes[x[si_nodes] <= 1e-12]
    drn = np.flatnonzero(x >= M.BOTTOM_DEPTH - 1e-9)

    for s in (body_nodes, gate, shield, src):
        fixed[s] = True
    val[fixed] = V_D
    fixed[drn] = True
    val[drn] = 0.0

    free = ~fixed
    idx = np.flatnonzero(free)

    # constraints on charged silicon only
    cons = np.flatnonzero(charged)
    nC = 6*len(cons)
    ids = tris[cons]
    rows = np.repeat(np.arange(nC), 6)
    cols = np.tile(ids, (1, 6)).ravel()
    vals = np.tile(B.ravel(), len(cons))
    C = sp.coo_matrix((vals, (rows, cols)), shape=(nC, len(nodes))).tocsr()
    C_F = C[:, idx].tocsc()
    cD = C[:, fixed] @ val[fixed]

    Kff = K[idx][:, idx].tocsc()
    g = -f[idx] - (K[idx][:, fixed] @ val[fixed])
    sK = float(Kff.diagonal().max())
    Kff_s = (Kff/sK).tocsc()
    g_s = g/sK

    print(f"[P2.3] dofs={len(nodes)} free={len(idx)} constraints={nC} "
          f"charged elements={len(cons)} body elements={int(body.sum())}")
    print(f"[P2.3] charge={A.RHO*area[charged].sum():.6f}")

    prob = osqp.OSQP()
    prob.setup(P=Kff_s, q=-g_s, A=C_F, l=-cD, u=np.full(len(cD), np.inf),
               verbose=False, eps_abs=1e-8, eps_rel=1e-8, max_iter=200000,
               polishing=True, scaling=10)
    res = prob.solve()
    print(f"[P2.3] OSQP {res.info.status} iters={res.info.iter} "
          f"obj={res.info.obj_val:.10e}")

    u = val.copy()
    u[idx] = res.x
    phi = V_D - u
    slack = C_F @ res.x + cD
    lam = -res.y
    r = Kff_s @ res.x - g_s
    stat = r - C_F.T @ lam

    print(f"[P2.3] min slack={slack.min():+.3e} active={int((slack < 1e-7).sum())}/{len(slack)}")
    print(f"[P2.3] lam min={lam.min():+.3e} max={lam.max():.3e}")
    print(f"[P2.3] stationarity={np.linalg.norm(stat)/max(np.linalg.norm(r),1e-30):.3e}")
    print(f"[P2.3] complementarity={np.max(np.abs(lam*slack)):.3e}")
    print(f"[P2.3] phi range=[{phi.min():.6f},{phi.max():.6f}] min u={u.min():.6e}")

    # neutral set: silicon elements with at least one active Bernstein row
    bs = u[tris[cons]] @ B.T
    act_el = np.zeros(len(tris), bool)
    act_el[cons] = (bs < 1e-7).any(axis=1)
    neutral_area = float(area[act_el].sum())
    si_area = float(area[si].sum())
    free_si = si_nodes[free[si_nodes]]
    imin = free_si[np.argmin(u[free_si])]
    print(f"[P2.3] neutral area={neutral_area:.9f} um^2 "
          f"({100*neutral_area/si_area:.4f}% of Si)")
    print(f"[P2.3] min u free Si={u[free_si].min():+.9f} V at "
          f"(x={nodes[imin,0]:.5f}, y={nodes[imin,1]:.5f})")

    np.savez_compressed("p2_p3_r3.npz", nodes=nodes, tris=tris, mat=mat,
                        phi=phi, u=u, slack=slack, lam=lam, active=act_el)
    json.dump(dict(status=str(res.info.status), iters=int(res.info.iter),
                   free=int(len(idx)), constraints=int(nC),
                   neutral_area=neutral_area, si_area=si_area,
                   min_u=float(u[free_si].min()),
                   active=int((slack < 1e-7).sum()),
                   stationarity=float(np.linalg.norm(stat)/max(np.linalg.norm(r), 1e-30)),
                   complementarity=float(np.max(np.abs(lam*slack)))),
              open("p2_p3_r3.json", "w"), indent=2)
    print("wrote p2_p3_r3.npz / p2_p3_r3.json")


if __name__ == "__main__":
    main()
