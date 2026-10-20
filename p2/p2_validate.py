#!/usr/bin/env python3
"""F7: hard full-set KKT gate on a polished P2 obstacle solution."""
from __future__ import annotations
import json
import numpy as np
import scipy.sparse as sp

import p2_audit as A
from p2_bernstein import B
from p2_nested_consistency import assemble_p2

M = A.M
V_D = A.V_DRAIN


def check(tag, fn):
    z = np.load(fn)
    nodes, tris, mat = z["nodes"], z["tris"], z["mat"]
    u, lam = z["u"], z["lam"]
    cxs = nodes[tris][:, :3, 0].mean(axis=1)
    si = mat == M.SI
    body = si & (cxs <= A.BODY_BOTTOM + 1e-12)
    charged = si & ~body
    K, f = assemble_p2(nodes, tris, mat, charged)
    x = nodes[:, 0]
    si_nodes = np.unique(tris[si])
    gate = np.unique(tris[mat == M.GATE]); shield = np.unique(tris[mat == M.SHIELD])
    src = si_nodes[x[si_nodes] <= 1e-12]; drn = np.flatnonzero(x >= M.BOTTOM_DEPTH - 1e-9)
    fixed = np.zeros(len(nodes), bool); val = np.zeros(len(nodes))
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
    uf = u[idx]
    slack = C_F @ uf + cD
    r = Kff @ uf - g
    stat = r - C_F.T @ lam
    comp = float(np.max(np.abs(lam*slack)))
    min_slack = float(slack.min())
    neg = int((slack < -1e-9).sum())
    lam_min = float(lam.min()); lam_max = float(lam.max())
    inact = np.abs(lam) < 1e-14
    inact_bad = int((np.abs(lam[inact]) > 0).sum())
    rel = float(np.linalg.norm(stat)/max(np.linalg.norm(r), 1e-30))
    overshoot = float(-u.min())
    si_set = set(np.unique(tris[si]).tolist())
    ox = [int(i) for i in np.unique(tris[mat == M.OX]) if int(i) not in si_set]
    ox_free = np.intersect1d(np.array(ox, dtype=int), idx)
    Ku = K @ u + f
    # Ku is full-length: index with GLOBAL node ids, not free-vector positions
    ox_res = float(np.abs(Ku[ox_free]).max()) if len(ox_free) else 0.0
    print(f"\n[{tag}] {fn}")
    print(f"  constraints={nC}  saved active={len(z['active'])}")
    print(f"  min slack={min_slack:+.3e}  violations<-1e-9={neg}")
    print(f"  lambda min/max={lam_min:+.3e}/{lam_max:+.3e}  inactive nonzero={inact_bad}")
    print(f"  stationarity rel={rel:.3e}  complementarity={comp:.3e}")
    print(f"  oxide residual={ox_res:.3e}  overshoot={overshoot:+.3e}")
    ok = (neg == 0 and lam_min >= -1e-9 and rel < 1e-8 and comp < 1e-10
          and ox_res < 1e-9 and overshoot < 1e-9)
    print(f"  GATE: {'PASS' if ok else 'FAIL'}")
    return dict(constraints=nC, active=int(len(z["active"])), min_slack=min_slack,
                violations=neg, lam_min=lam_min, lam_max=lam_max, inact_nonzero=inact_bad,
                stationarity=rel, complementarity=comp, oxide=ox_res,
                overshoot=overshoot, gate=bool(ok))


if __name__ == "__main__":
    out = {}
    for tag, fn in (("r3", "p2_p3_r3_polished.npz"),
                    ("r4", "p2_p3_r4_polished.npz"),
                    ("mid", "p2_p3_rmid_polished.npz")):
        out[tag] = check(tag, fn)
    json.dump(out, open("p2_validate.json", "w"), indent=2)
    print("\nwrote p2_validate.json")
