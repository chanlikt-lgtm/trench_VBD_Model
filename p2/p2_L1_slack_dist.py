#!/usr/bin/env python3
"""Size the L1 active set: slack and multiplier distributions at the OSQP point."""
from __future__ import annotations
import numpy as np
import scipy.sparse as sp

import p2_audit as A
from p2_bernstein import B
from p2_nested_consistency import assemble_p2

M = A.M


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
    val[fixed] = 34.274208536
    fixed[drn] = True; val[drn] = 0.0
    idx = np.flatnonzero(~fixed)
    cons = np.flatnonzero(charged); nC = 6*len(cons)
    C = sp.coo_matrix((np.tile(B.ravel(), len(cons)),
                       (np.repeat(np.arange(nC), 6), np.tile(tris[cons], (1, 6)).ravel())),
                      shape=(nC, len(nodes))).tocsr()
    C_F = C[:, idx].tocsc()
    cD = C[:, fixed] @ val[fixed]
    chk = np.load("p2_L1_osqp_checkpoint.npz")
    xf = chk["x"]; y = chk["y"]; sK = float(chk["sK"])
    lam = -sK*y
    slack = C_F @ xf + cD
    print(f"min slack={slack.min():+.3e}  max slack={slack.max():.3e}")
    print("slack counts:")
    for th in (1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7):
        print(f"  <{th:g}: {int((slack < th).sum())}")
    print("lambda counts:")
    for th in (1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4):
        print(f"  >{th:g}: {int((lam > th).sum())}")
    for sth in (1e-10, 1e-9, 1e-8):
        for lth in (1e-10, 1e-8, 1e-6):
            print(f"  (slack<{sth:g} & lam>{lth:g}): "
                  f"{int(((slack < sth) & (lam > lth)).sum())}")


if __name__ == "__main__":
    main()
