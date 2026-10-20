#!/usr/bin/env python3
"""F3/F9 cleanup: edge-connected neutral topology and threshold-resolved cuts.

- Components of the fully neutral set using shared EDGES (two common vertices),
  compared with the previous shared-node criterion.
- Front intersections at x = 2.35 and 2.41 for eps in {1e-8, 1e-7, 1e-6, 1e-5},
  merging adjacent element segments and separating short artifacts.
"""
from __future__ import annotations
import json
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

import p2_audit as A
from p2_exact_metrics import cut_intervals

M = A.M


def components(nodes, tris, mat, u, mode):
    ue = u[tris]
    cxs = nodes[tris][:, :3, 0].mean(axis=1)
    charged = (mat == M.SI) & (cxs > A.BODY_BOTTOM + 1e-12)
    fully = charged & (np.abs(ue) < 1e-8).all(axis=1)
    els = np.flatnonzero(fully)
    local = {int(e): i for i, e in enumerate(els)}
    buckets = {}
    for e in els:
        t = tris[e, :3]
        keys = []
        if mode == "edge":
            keys = [(min(t[0], t[1]), max(t[0], t[1])),
                    (min(t[1], t[2]), max(t[1], t[2])),
                    (min(t[2], t[0]), max(t[2], t[0]))]
        else:
            keys = [int(v) for v in t]
        for k in keys:
            buckets.setdefault(k, []).append(int(e))
    rows, cols = [], []
    for k, es in buckets.items():
        for a in es:
            for b in es:
                if a != b:
                    rows.append(local[a]); cols.append(local[b])
    Adj = coo_matrix((np.ones(len(rows)), (rows, cols)),
                     shape=(len(els), len(els))).tocsr()
    ncomp, lab = connected_components(Adj, directed=False) if len(els) else (0, np.array([]))
    area = np.abs(A.element_geom(nodes, tris)[2])
    drain = nodes[:, 0] >= M.BOTTOM_DEPTH - 1e-9
    comps = []
    for c in range(ncomp):
        m = lab == c
        es = els[m]
        pts = nodes[np.unique(tris[es])]
        comps.append(dict(n=int(m.sum()), area=float(area[es].sum()),
                          x=[float(pts[:, 0].min()), float(pts[:, 0].max())],
                          y=[float(pts[:, 1].min()), float(pts[:, 1].max())],
                          drain=bool(drain[np.unique(tris[es])].any())))
    return dict(total=int(fully.sum()), components=ncomp, comps=comps)


def clean(iv, merge_tol=1e-3, min_len=1e-3):
    if not iv:
        return [], []
    iv = sorted(iv)
    merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1] + merge_tol:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    main = [(a, b) for a, b in merged if b-a >= min_len]
    short = [(a, b) for a, b in merged if b-a < min_len]
    return main, short


def main():
    out = {}
    for tag, fn in (("r3", "p2_p3_r3_polished.npz"),
                    ("r4", "p2_p3_r4_polished.npz"),
                    ("mid", "p2_p3_rmid_polished.npz")):
        z = np.load(fn)
        nodes, tris, mat, u = z["nodes"], z["tris"], z["mat"], z["u"]
        ce = components(nodes, tris, mat, u, "edge")
        cn = components(nodes, tris, mat, u, "node")
        print(f"\n[{tag}] fully neutral {ce['total']} elements")
        print(f"  edge-connected components: {ce['components']}")
        for c in ce["comps"]:
            print(f"    n={c['n']:4d} area={c['area']:.9f} x={c['x']} y={c['y']} "
                  f"drain={c['drain']}")
        print(f"  node-connected components (old): {cn['components']}")

        cuts = {}
        for x0 in (2.35, 2.41):
            cuts[str(x0)] = {}
            for eps in (1e-8, 1e-7, 1e-6, 1e-5):
                iv = cut_intervals(nodes, tris, mat, u, x0, eps, True)
                main_, short_ = clean(iv)
                cuts[str(x0)][str(eps)] = dict(main=main_, short_count=len(short_))
                print(f"  x={x0:.2f} eps={eps:g}: main={[(round(a,4),round(b,4)) for a,b in main_]}"
                      f"  short_artifacts={len(short_)}")
        out[tag] = dict(edge=ce, node=cn, cuts=cuts)
    json.dump(out, open("p2_topology_clean.json", "w"), indent=2)
    print("\nwrote p2_topology_clean.json")


if __name__ == "__main__":
    main()
