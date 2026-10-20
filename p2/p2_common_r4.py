#!/usr/bin/env python3
"""Common-criterion neutral area and front cuts for the r4 obstacle solution."""
from __future__ import annotations
import json
import sys
import numpy as np
import p2_audit as A
from p2_common_zero import sample_points

M = A.M


def main():
    FN = sys.argv[1] if len(sys.argv) > 1 else "p2_p3_r4_polished.npz"
    OUT = sys.argv[2] if len(sys.argv) > 2 else "p2_common_zero_r4.json"
    z = np.load(FN)
    nodes, tris, mat, u = z["nodes"], z["tris"], z["mat"], z["u"]
    area = np.abs(A.element_geom(nodes, tris)[2])
    cxs = nodes[tris][:, :3, 0].mean(axis=1)
    charged = (mat == M.SI) & (cxs > A.BODY_BOTTOM + 1e-12)
    els = np.flatnonzero(charged)

    L1, L2 = sample_points(0.05)
    N2 = np.array([A.shape_vals(a, b) for a, b in zip(L1, L2)])
    print(f"[r4 common] points={len(L1)} charged={len(els)}")
    out = {}
    for eps in (1e-8, 1e-6, 1e-4):
        a = 0.0
        for s in range(0, len(els), 2000):
            b = els[s:s+2000]
            U = u[tris[b]] @ N2.T
            a += float(((U <= eps).mean(axis=1)*area[b]).sum())
        out[str(eps)] = a
        print(f"  eps={eps:g}: neutral area={a:.9f}")

    # front cuts
    Aa = nodes[tris[els, 0]]; Bb = nodes[tris[els, 1]]; Cc = nodes[tris[els, 2]]
    det = ((Bb[:, 0]-Aa[:, 0])*(Cc[:, 1]-Aa[:, 1]) -
           (Cc[:, 0]-Aa[:, 0])*(Bb[:, 1]-Aa[:, 1]))

    def cut(x0, ys):
        vals = []
        for y in ys:
            l1 = ((Bb[:, 1]-Cc[:, 1])*(x0-Cc[:, 0]) +
                  (Cc[:, 0]-Bb[:, 0])*(y-Cc[:, 1]))/det
            l2 = ((Cc[:, 1]-Aa[:, 1])*(x0-Cc[:, 0]) +
                  (Aa[:, 0]-Cc[:, 0])*(y-Cc[:, 1]))/det
            l3 = 1-l1-l2
            m = (l1 >= -1e-9) & (l2 >= -1e-9) & (l3 >= -1e-9)
            if not m.any():
                vals.append(np.nan); continue
            j = int(np.flatnonzero(m)[0]); e = int(els[j])
            vals.append(float(np.dot(u[tris[e]], A.shape_vals(l1[j], l2[j]))))
        return np.array(vals)

    ys = np.arange(-0.6, 0.6001, 0.002)
    cuts = {}
    print("\n[r4 front cuts] eps=1e-8")
    for x0 in (2.35, 2.38, 2.39, 2.40, 2.41):
        v = cut(x0, ys)
        neu = (~np.isnan(v)) & (np.abs(v) <= 1e-8)
        iv = []
        if neu.any():
            idx = np.flatnonzero(neu); s = idx[0]; p = idx[0]
            for k in idx[1:]:
                if k != p+1:
                    iv.append((float(ys[s]), float(ys[p]))); s = k
                p = k
            iv.append((float(ys[s]), float(ys[p])))
        cuts[x0] = iv
        print(f"  x={x0:.2f}: {[(round(a,3), round(b,3)) for a,b in iv]}")

    json.dump(dict(areas=out, cuts={str(k): v for k, v in cuts.items()}),
              open(OUT, "w"), indent=2)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
