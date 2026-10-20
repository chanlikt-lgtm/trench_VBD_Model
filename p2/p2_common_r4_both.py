#!/usr/bin/env python3
"""Common-criterion P1/P2 comparison at r4."""
from __future__ import annotations
import json
import numpy as np
import p2_audit as A
from p2_common_zero import sample_points

M = A.M


def main():
    p1 = np.load("p1_on_p2mesh_obstacle_r4.npz")
    p2 = np.load("p2_p3_r4_polished.npz")
    nodes6, tris6, mat = p2["nodes"], p2["tris"], p2["mat"]
    u1, u2 = p1["u"], p2["u"]
    tv = tris6[:, :3]
    used = np.unique(tv.ravel())
    o2n = np.full(len(nodes6), -1, dtype=int)
    o2n[used] = np.arange(len(used))
    tris1 = o2n[tv]
    area = np.abs(A.element_geom(nodes6, tris6)[2])
    cxs = nodes6[tris6][:, :3, 0].mean(axis=1)
    charged = (mat == M.SI) & (cxs > A.BODY_BOTTOM + 1e-12)
    els = np.flatnonzero(charged)

    L1, L2 = sample_points(0.05)
    N1 = np.column_stack([L1, L2, 1-L1-L2])
    N2 = np.array([A.shape_vals(a, b) for a, b in zip(L1, L2)])
    print(f"[r4 common] points={len(L1)} charged={len(els)}")
    out = {}
    for eps in (1e-8, 1e-6, 1e-4):
        a1 = a2 = 0.0
        for s in range(0, len(els), 2000):
            b = els[s:s+2000]
            U1 = u1[tris1[b]] @ N1.T
            U2 = u2[tris6[b]] @ N2.T
            a1 += float(((U1 <= eps).mean(axis=1)*area[b]).sum())
            a2 += float(((U2 <= eps).mean(axis=1)*area[b]).sum())
        out[str(eps)] = [a1, a2]
        print(f"  eps={eps:g}: P1={a1:.9f}  P2={a2:.9f}  "
              f"ratio={a2/a1:.4f}  diff={100*(a2-a1)/a1:+.2f}%")

    def cut(field, nodes_src, tris_src, x0, ys, quad):
        cidx = els
        Aa = nodes_src[tris_src[cidx, 0]]; Bb = nodes_src[tris_src[cidx, 1]]
        Cc = nodes_src[tris_src[cidx, 2]]
        det = ((Bb[:, 0]-Aa[:, 0])*(Cc[:, 1]-Aa[:, 1]) -
               (Cc[:, 0]-Aa[:, 0])*(Bb[:, 1]-Aa[:, 1]))
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
            j = int(np.flatnonzero(m)[0]); e = int(cidx[j])
            if quad:
                vals.append(float(np.dot(field[tris6[e]], A.shape_vals(l1[j], l2[j]))))
            else:
                vals.append(float(l1[j]*field[tris1[e, 0]] +
                                  l2[j]*field[tris1[e, 1]] +
                                  (1-l1[j]-l2[j])*field[tris1[e, 2]]))
        return np.array(vals)

    ys = np.arange(-0.6, 0.6001, 0.002)
    cuts = {}
    print("\n[r4 front cuts] eps=1e-8")
    for x0 in (2.35, 2.38, 2.39, 2.40, 2.41):
        v1 = cut(u1, p1["nodes"], tris1, x0, ys, False)
        v2 = cut(u2, nodes6, tris6, x0, ys, True)
        def iv(v):
            neu = (~np.isnan(v)) & (np.abs(v) <= 1e-8)
            if not neu.any():
                return []
            idx = np.flatnonzero(neu); res = []; s = idx[0]; p = idx[0]
            for k in idx[1:]:
                if k != p+1:
                    res.append((float(ys[s]), float(ys[p]))); s = k
                p = k
            res.append((float(ys[s]), float(ys[p])))
            return res
        cuts[x0] = dict(p1=iv(v1), p2=iv(v2))
        print(f"  x={x0:.2f}: P1={[(round(a,3), round(b,3)) for a,b in cuts[x0]['p1']]}")
        print(f"           P2={[(round(a,3), round(b,3)) for a,b in cuts[x0]['p2']]}")

    json.dump(dict(areas=out, cuts={str(k): v for k, v in cuts.items()}),
              open("p2_common_zero_r4.json", "w"), indent=2)
    print("\nwrote p2_common_zero_r4.json")


if __name__ == "__main__":
    main()
