#!/usr/bin/env python3
"""Step B: full-set F7 gate on the saved L1 OSQP solution."""
from __future__ import annotations
import json
import numpy as np
import scipy.sparse as sp

import p2_audit as A
from p2_bernstein import B
from p2_nested_consistency import assemble_p2

M = A.M
V_D = A.V_DRAIN


def main():
    z = np.load("p2_L1_p2mesh.npz")
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
    K, f = assemble_p2(nodes, tris, mat, charged)

    fixed = np.zeros(len(nodes), bool); val = np.zeros(len(nodes))
    si_nodes = np.unique(tris[si])
    gate = np.unique(tris[mat == M.GATE]); shield = np.unique(tris[mat == M.SHIELD])
    src = si_nodes[x[si_nodes] <= 1e-12]; drn = np.flatnonzero(x >= M.BOTTOM_DEPTH - 1e-9)
    for s in (np.unique(tris[body]), gate, shield, src):
        fixed[s] = True
    val[fixed] = V_D
    fixed[drn] = True; val[drn] = 0.0
    idx = np.flatnonzero(~fixed)

    cons = np.flatnonzero(charged); nC = 6*len(cons)
    C = sp.coo_matrix((np.tile(B.ravel(), len(cons)),
                       (np.repeat(np.arange(nC), 6), np.tile(tris[cons], (1, 6)).ravel())),
                      shape=(nC, len(nodes))).tocsr()
    C_F = C[:, idx].tocsc()
    cD = C[:, fixed] @ val[fixed]
    Kff = K[idx][:, idx].tocsc()
    g = -f[idx] - (K[idx][:, fixed] @ val[fixed])

    chk = np.load("p2_L1_osqp_checkpoint.npz")
    xf = chk["x"]; y = chk["y"]; sK = float(chk["sK"])
    lam = -sK*y
    slack = C_F @ xf + cD
    r = Kff @ xf - g
    stat = r - C_F.T @ lam
    comp = float(np.max(np.abs(lam*slack)))
    rel = float(np.linalg.norm(stat)/max(np.linalg.norm(r), 1e-30))

    u = val.copy(); u[idx] = xf
    overshoot = float(-u.min())
    si_set = set(np.unique(tris[si]).tolist())
    ox = [int(i) for i in np.unique(tris[mat == M.OX]) if int(i) not in si_set]
    ox_free = np.intersect1d(np.array(ox, dtype=int), idx)
    Ku = K @ u + f
    ox_res = float(np.abs(Ku[ox_free]).max()) if len(ox_free) else 0.0

    print(f"[F7-OSQP] constraints={nC}  free={len(idx)}")
    print(f"  min slack={slack.min():+.3e}  violations<-1e-9={int((slack<-1e-9).sum())}")
    print(f"  lambda min/max={lam.min():+.3e}/{lam.max():+.3e}")
    print(f"  stationarity rel={rel:.3e}")
    print(f"  complementarity={comp:.3e}")
    print(f"  oxide residual={ox_res:.3e}")
    print(f"  overshoot={overshoot:+.3e}")
    ok = (int((slack < -1e-9).sum()) == 0 and lam.min() >= -1e-9 and
          rel < 1e-8 and comp < 1e-10 and ox_res < 1e-9 and overshoot < 1e-9)
    print(f"  GATE: {'PASS' if ok else 'FAIL'}")
    json.dump(dict(constraints=nC, free=int(len(idx)), min_slack=float(slack.min()),
                   violations=int((slack < -1e-9).sum()), lam_min=float(lam.min()),
                   lam_max=float(lam.max()), stationarity=rel, complementarity=comp,
                   oxide=ox_res, overshoot=overshoot, gate=bool(ok)),
              open("p2_L1_f7_osqp.json", "w"), indent=2)
    print("wrote p2_L1_f7_osqp.json")


if __name__ == "__main__":
    main()
