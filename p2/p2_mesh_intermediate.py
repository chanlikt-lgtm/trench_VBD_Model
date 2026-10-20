#!/usr/bin/env python3
"""Mesh-only step for the intermediate P2 resolution (h=0.0085, h_tip=0.0034).

Builds the conforming mesh, verifies quality, and reports the exact DOF and
constraint counts plus a memory estimate. No assembly, no OSQP.
"""
from __future__ import annotations
import json
import numpy as np
import p2_build as B
import p2_audit as A

M = A.M
H, H_TIP = 0.0085, 0.0034


def main():
    nodes, tris, mat, counts = B.build_p2(H, H_TIP, "p2_rmid")
    n = len(nodes)
    used = np.unique(tris.ravel())
    n_used = len(used)
    p = nodes[tris][:, :3]
    det = ((p[:, 1, 0]-p[:, 0, 0])*(p[:, 2, 1]-p[:, 0, 1]) -
           (p[:, 2, 0]-p[:, 0, 0])*(p[:, 1, 1]-p[:, 0, 1]))
    area = 0.5*np.abs(det)

    cxs = nodes[tris][:, :3, 0].mean(axis=1)
    si = mat == M.SI
    body = si & (cxs <= A.BODY_BOTTOM + 1e-12)
    charged = si & ~body
    x = nodes[:, 0]
    si_nodes = np.unique(tris[si])
    gate = np.unique(tris[mat == M.GATE])
    shield = np.unique(tris[mat == M.SHIELD])
    src = si_nodes[x[si_nodes] <= 1e-12]
    drn = np.flatnonzero(x >= M.BOTTOM_DEPTH - 1e-9)
    fixed = np.zeros(n, bool)
    for s in (np.unique(tris[body]), gate, shield, src):
        fixed[s] = True
    fixed[drn] = True
    free = int((~fixed).sum())
    ncons = 6*int(charged.sum())

    # quality
    edges = np.stack([np.linalg.norm(p[:, 1]-p[:, 0], axis=1),
                      np.linalg.norm(p[:, 2]-p[:, 1], axis=1),
                      np.linalg.norm(p[:, 0]-p[:, 2], axis=1)], 1)
    emax = edges.max(axis=1)
    hmin = 2*area/emax
    aspect = emax/np.maximum(hmin, 1e-30)
    ang = []
    for k in range(3):
        v1 = p[:, (k+1) % 3]-p[:, k]
        v2 = p[:, (k+2) % 3]-p[:, k]
        cs = np.sum(v1*v2, 1)/(np.linalg.norm(v1, axis=1)*np.linalg.norm(v2, axis=1)+1e-30)
        ang.append(np.degrees(np.arccos(np.clip(cs, -1, 1))))
    ang = np.concatenate(ang)

    # memory estimate
    c_nnz = 12*int(charged.sum())
    k_nnz_lo = 18*len(tris)
    k_nnz_hi = 24*len(tris)
    print(f"[mesh] nodes={n} (used {n_used}) elements={len(tris)}")
    print(f"  Si={int(si.sum())} OX={int((mat==M.OX).sum())} "
          f"Gate={int((mat==M.GATE).sum())} Shield={int((mat==M.SHIELD).sum())}")
    print(f"  charged={int(charged.sum())} body={int(body.sum())}")
    print(f"  fixed={int(fixed.sum())} free={free} constraints={ncons}")
    print(f"  area total={area.sum():.6f} Si={area[si].sum():.6f} "
          f"charged={area[charged].sum():.6f}")
    print(f"  quality: min area={area.min():.3e} max aspect={aspect.max():.2f} "
          f"min angle={ang.min():.2f} deg")
    print(f"  estimate: K nnz {k_nnz_lo/1e6:.2f}-{k_nnz_hi/1e6:.2f}M  "
          f"C nnz {c_nnz/1e6:.3f}M")
    print(f"  r4 reference: 296,014 DOFs, 240,606 free, 583,674 constraints")

    out = dict(h=H, h_tip=H_TIP, nodes=int(n), nodes_used=int(n_used),
               elements=int(len(tris)), si=int(si.sum()), ox=int((mat == M.OX).sum()),
               gate=int((mat == M.GATE).sum()), shield=int((mat == M.SHIELD).sum()),
               charged=int(charged.sum()), body=int(body.sum()),
               fixed=int(fixed.sum()), free=free, constraints=ncons,
               area_total=float(area.sum()), area_si=float(area[si].sum()),
               area_charged=float(area[charged].sum()),
               min_area=float(area.min()), max_aspect=float(aspect.max()),
               min_angle=float(ang.min()),
               K_nnz_lo=int(k_nnz_lo), K_nnz_hi=int(k_nnz_hi), C_nnz=int(c_nnz))
    json.dump(out, open("p2_rmid_report.json", "w"), indent=2)
    print("wrote p2_rmid_report.json")


if __name__ == "__main__":
    main()
