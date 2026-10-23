#!/usr/bin/env python3
"""M2 Step A: solve OSQP and checkpoint immediately (no polish)."""
from __future__ import annotations
import json
import numpy as np
import scipy.sparse as sp
import osqp

import p2_audit as A
from p2_bernstein import B
from p2_nested_consistency import assemble_p2
from p2_exact_donor import exact_load

M = A.M
V_D = A.V_DRAIN


def main():
    z = np.load("p2_B2_p2mesh.npz")
    nodes, tris, mat = z["nodes"], z["tris"], z["mat"]
    used = np.unique(tris.ravel())
    o2n = np.full(len(nodes), -1, dtype=int)
    o2n[used] = np.arange(len(used))
    nodes = nodes[used]; tris = o2n[tris]
    x = nodes[:, 0]

    cxs = nodes[tris][:, :3, 0].mean(axis=1)
    si = mat == M.SI
    body = si & (cxs <= A.BODY_BOTTOM + 1e-12)
    charged = si & ~body
    K, _f0 = assemble_p2(nodes, tris, mat, np.zeros(len(tris), bool))
    f, _ftot = exact_load(nodes, tris, mat)

    fixed = np.zeros(len(nodes), bool)
    val = np.zeros(len(nodes))
    si_nodes = np.unique(tris[si])
    gate = np.unique(tris[mat == M.GATE]); shield = np.unique(tris[mat == M.SHIELD])
    src = si_nodes[x[si_nodes] <= 1e-12]; drn = np.flatnonzero(x >= M.BOTTOM_DEPTH - 1e-9)
    for s in (np.unique(tris[body]), gate, shield, src):
        fixed[s] = True
    val[fixed] = V_D
    fixed[drn] = True; val[drn] = 0.0
    idx = np.flatnonzero(~fixed)

    cons = np.flatnonzero(charged)
    nC = 6*len(cons)
    C = sp.coo_matrix((np.tile(B.ravel(), len(cons)),
                       (np.repeat(np.arange(nC), 6), np.tile(tris[cons], (1, 6)).ravel())),
                      shape=(nC, len(nodes))).tocsr()
    C_F = C[:, idx].tocsc()
    cD = C[:, fixed] @ val[fixed]
    Kff = K[idx][:, idx].tocsc()
    g = -f[idx] - (K[idx][:, fixed] @ val[fixed])
    sK = float(Kff.diagonal().max())
    Kff_s = (Kff/sK).tocsc(); g_s = g/sK

    print(f"[B2X-A] dofs={len(nodes)} free={len(idx)} constraints={nC} "
          f"charged={len(cons)} body={int(body.sum())}", flush=True)
    prob = osqp.OSQP()
    prob.setup(P=Kff_s, q=-g_s, A=C_F, l=-cD, u=np.full(len(cD), np.inf),
               verbose=False, eps_abs=1e-8, eps_rel=1e-8, max_iter=30000,
               polishing=False, scaling=10)
    res = prob.solve()
    pri = float(getattr(res.info, "prim_res", getattr(res.info, "pri_res", float("nan"))))
    dua = float(getattr(res.info, "dual_res", getattr(res.info, "dua_res", float("nan"))))
    print(f"[B2X-A] status={res.info.status} iters={res.info.iter} "
          f"obj={res.info.obj_val:.6e} pri={pri:.3e} dua={dua:.3e}", flush=True)

    np.savez_compressed("p2_B2X_osqp_checkpoint.npz",
                        x=res.x, y=res.y, idx=idx, fixed=fixed, val=val,
                        cons=cons, status=str(res.info.status),
                        iters=int(res.info.iter), obj=float(res.info.obj_val),
                        pri_res=pri, dua_res=dua,
                        sK=sK)
    json.dump(dict(status=str(res.info.status), iters=int(res.info.iter),
                   obj=float(res.info.obj_val), pri_res=pri,
                   dua_res=dua, sK=sK,
                   dofs=int(len(nodes)), free=int(len(idx)), constraints=int(nC),
                   charged=int(len(cons)), body=int(body.sum())),
              open("p2_B2X_osqp.json", "w"), indent=2)
    # verify reload
    chk = np.load("p2_B2X_osqp_checkpoint.npz")
    print(f"[B2X-A] checkpoint reload: x{chk['x'].shape} y{chk['y'].shape} "
          f"idx{chk['idx'].shape} status={chk['status']}", flush=True)
    print("wrote p2_B2X_osqp_checkpoint.npz / p2_B2X_osqp.json")


if __name__ == "__main__":
    main()


