#!/usr/bin/env python3
"""CP-M2-PRE1..PRE4 handoff diagnostics for the frozen S0-L0/L1 pair.

Read-only: no assembly, no solve, no area estimation.
"""
from __future__ import annotations
import json
import numpy as np
from scipy.spatial import cKDTree

import p2_audit as A

M = A.M


def load(tag):
    z = np.load(f"../p2_3_spatial_S0/levels/S0_{tag}_p2mesh.npz")
    return z["nodes"], z["tris"], z["mat"]


def boundary_nodes(nodes, tris, mat, code):
    """Nodes on the boundary of a material region."""
    sel = tris[mat == code][:, :3]
    edges = np.sort(np.vstack([sel[:, [0, 1]], sel[:, [1, 2]], sel[:, [2, 0]]]), axis=1)
    uniq, cnt = np.unique(edges, axis=0, return_counts=True)
    b = uniq[cnt == 1]
    return np.unique(b.ravel())


def midside_dev(nodes, tris):
    P = nodes[tris]
    v0, v1, v2 = P[:, 0], P[:, 1], P[:, 2]
    m01, m12, m20 = P[:, 3], P[:, 4], P[:, 5]
    d = np.max([np.abs(m01 - 0.5*(v0+v1)).max(),
                np.abs(m12 - 0.5*(v1+v2)).max(),
                np.abs(m20 - 0.5*(v2+v0)).max()])
    return float(d)


def orientation(nodes, tris, mat):
    P = nodes[tris][:, :3]
    det = ((P[:, 1, 0]-P[:, 0, 0])*(P[:, 2, 1]-P[:, 0, 1]) -
           (P[:, 2, 0]-P[:, 0, 0])*(P[:, 1, 1]-P[:, 0, 1]))
    out = {}
    for code in np.unique(mat):
        m = mat == code
        out[str(int(code))] = dict(n=int(m.sum()),
                                   pos=int((det[m] > 0).sum()),
                                   neg=int((det[m] < 0).sum()),
                                   zero=int((det[m] == 0).sum()))
    return out


def quality(nodes, tris):
    P = nodes[tris][:, :3]
    det = ((P[:, 1, 0]-P[:, 0, 0])*(P[:, 2, 1]-P[:, 0, 1]) -
           (P[:, 2, 0]-P[:, 0, 0])*(P[:, 1, 1]-P[:, 0, 1]))
    area = 0.5*np.abs(det)
    edges = np.stack([np.linalg.norm(P[:, 1]-P[:, 0], axis=1),
                      np.linalg.norm(P[:, 2]-P[:, 1], axis=1),
                      np.linalg.norm(P[:, 0]-P[:, 2], axis=1)], 1)
    emax = edges.max(1); hmin = 2*area/np.maximum(emax, 1e-30)
    aspect = emax/np.maximum(hmin, 1e-30)
    ang = []
    for k in range(3):
        v1 = P[:, (k+1) % 3]-P[:, k]; v2 = P[:, (k+2) % 3]-P[:, k]
        cs = np.sum(v1*v2, 1)/(np.linalg.norm(v1, axis=1)*np.linalg.norm(v2, axis=1)+1e-30)
        ang.append(np.degrees(np.arccos(np.clip(cs, -1, 1))))
    ang = np.concatenate(ang)
    code = None
    return dict(min_angle=float(ang.min()),
                p01_angle=float(np.percentile(ang, 1)),
                p05_angle=float(np.percentile(ang, 5)),
                aspect_max=float(aspect.max()),
                aspect_p99=float(np.percentile(aspect, 99)),
                n_aspect_gt3=int((aspect > 3).sum()),
                n_aspect_gt5=int((aspect > 5).sum()),
                n_inverted=int((det <= 0).sum()),
                min_area=float(area.min()))


def donor_identity(nodes, tris, mat):
    from p2_exact_donor import exact_load
    from p2_charge_mask_audit import clip_area
    f, tot = exact_load(nodes, tris, mat)
    P = nodes[tris][:, :3]
    area = np.abs(A.element_geom(nodes, tris)[2])
    si = mat == M.SI
    xmin = P[:, :, 0].min(1); xmax = P[:, :, 0].max(1)
    exact = area[si & (xmin > A.BODY_BOTTOM)].sum()
    for k in np.flatnonzero(si & (xmin <= A.BODY_BOTTOM) & (xmax > A.BODY_BOTTOM)):
        exact += clip_area([tuple(q) for q in P[k]], A.BODY_BOTTOM)
    return float(tot), float(exact), float(tot-exact)


def main():
    out = {}
    data = {t: load(t) for t in ("L0", "L1")}

    for t, (nodes, tris, mat) in data.items():
        P = nodes[tris][:, :3]
        cx = P[:, :, 0].mean(1)
        charged = (mat == M.SI) & (cx > A.BODY_BOTTOM + 1e-12)
        used = np.unique(tris.ravel())
        tot, exact, d = donor_identity(nodes, tris, mat)
        out[t] = dict(
            nodes=int(len(nodes)),
            elements=int(len(tris)),
            p2_dofs_raw=int(len(nodes)),
            charged_elements=int(charged.sum()),
            bernstein_rows=int(6*charged.sum()),
            orphan_nodes=int(len(nodes)-len(used)),
            midside_max_dev=midside_dev(nodes, tris),
            quality=quality(nodes, tris),
            p2_dofs_compacted=int(len(used)),
            orientation=orientation(nodes, tris, mat),
            donor_sum=tot, donor_exact=exact, donor_delta=d,
        )
        print(f"[{t}] midside dev {out[t]['midside_max_dev']:.2e}  "
              f"min angle {out[t]['quality']['min_angle']:.2f}  "
              f"aspect max {out[t]['quality']['aspect_max']:.3f}  "
              f"orphans {out[t]['orphan_nodes']}")
        print(f"     dofs raw {out[t]['p2_dofs_raw']} compacted {out[t]['p2_dofs_compacted']}"
              f"  charged {out[t]['charged_elements']}  rows {out[t]['bernstein_rows']}"
              f"  donor delta {d:+.2e}")
        for k, v in out[t]["orientation"].items():
            print(f"     orientation mat {k}: n={v['n']} pos={v['pos']} neg={v['neg']} zero={v['zero']}")

    # PRE1: boundary drift of L1 nodes against L0 boundary segments
    n0, t0, m0 = data["L0"]; n1, t1, m1 = data["L1"]
    bn0 = boundary_nodes(n0, t0, m0, M.SI)
    bn1 = boundary_nodes(n1, t1, m1, M.SI)
    sel = t0[m0 == M.SI][:, :3]
    edges = np.sort(np.vstack([sel[:, [0, 1]], sel[:, [1, 2]], sel[:, [2, 0]]]), axis=1)
    uniq, cnt = np.unique(edges, axis=0, return_counts=True)
    bnd = uniq[cnt == 1]
    Aseg = n0[bnd[:, 0]]; Bseg = n0[bnd[:, 1]]
    pts = n1[bn1]
    best = np.full(len(pts), np.inf)
    for i in range(len(Aseg)):
        a = Aseg[i]; ab = Bseg[i]-a; L2 = float(ab@ab)
        t = np.clip(((pts-a) @ ab)/L2, 0.0, 1.0)
        proj = a + t[:, None]*ab
        d = np.linalg.norm(pts-proj, axis=1)
        best = np.minimum(best, d)
    d = best
    out["PRE1"] = dict(l0_boundary_nodes=int(len(bn0)), l1_boundary_nodes=int(len(bn1)),
                       max_l1_node_to_l0_segment=float(d.max()),
                       p99=float(np.percentile(d, 99)))
    print(f"[PRE1] L1 boundary nodes vs L0 segments: max {d.max():.3e} um "
          f"(p99 {np.percentile(d,99):.3e})")

    json.dump(out, open("../p2_3_spatial_S0/CP-M2-PRE.json", "w"), indent=2)
    print("wrote CP-M2-PRE.json")


if __name__ == "__main__":
    main()
