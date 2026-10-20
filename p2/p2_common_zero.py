#!/usr/bin/env python3
"""Common P1/P2 neutral-area definition.

Both fields live on the same vertex triangles. The neutral region is the
sub-level set {u <= eps} with eps = 1e-8 V, and its area is computed by dense
barycentric sampling inside each charged silicon element, using each field's
own interpolation (linear for P1, quadratic for P2). This replaces the native
active-node / fully-neutral-element counts with one geometric criterion.
"""
from __future__ import annotations
import json
import numpy as np
import p2_audit as A

M = A.M
V_D = A.V_DRAIN
EPS_LIST = (1e-8, 1e-6, 1e-4)


def sample_points(step=0.05):
    L1s, L2s = [], []
    n = int(round(1.0/step))
    for i in range(n+1):
        for j in range(n+1-i):
            L1s.append(i*step); L2s.append(j*step)
    L1 = np.array(L1s); L2 = np.array(L2s)
    m = (L1+L2) <= 1.0 + 1e-12
    return L1[m], L2[m]


def main():
    p1 = np.load("p1_on_p2mesh_obstacle.npz")
    p2 = np.load("p2_p3_r3_polished.npz")
    nodes6, tris6, mat = p2["nodes"], p2["tris"], p2["mat"]
    u1 = p1["u"]; u2 = p2["u"]
    # P1 solution is on the compacted vertex mesh: map original P2 vertices
    tv = tris6[:, :3]
    used = np.unique(tv.ravel())
    o2n = np.full(len(nodes6), -1, dtype=int)
    o2n[used] = np.arange(len(used))
    tris1 = o2n[tv]
    area = np.abs(A.element_geom(nodes6, tris6)[2])
    cxs = nodes6[tris6][:, :3, 0].mean(axis=1)
    si = mat == M.SI
    charged = si & (cxs > A.BODY_BOTTOM + 1e-12)
    els = np.flatnonzero(charged)

    L1, L2 = sample_points(0.05)
    N1 = np.column_stack([L1, L2, 1-L1-L2])
    N2 = np.array([A.shape_vals(a, b) for a, b in zip(L1, L2)])
    print(f"[common] sample points per element={len(L1)}  charged elements={len(els)}")

    out = {}
    for eps in EPS_LIST:
        a1 = a2 = 0.0
        for s in range(0, len(els), 2000):
            b = els[s:s+2000]
            U1 = u1[tris1[b]] @ N1.T
            U2 = u2[tris6[b]] @ N2.T
            a1 += float(((U1 <= eps).mean(axis=1) * area[b]).sum())
            a2 += float(((U2 <= eps).mean(axis=1) * area[b]).sum())
        out[eps] = (a1, a2)
        print(f"  eps={eps:g}: P1 neutral area={a1:.9f}  P2={a2:.9f}  "
              f"ratio P2/P1={a2/a1:.4f}  diff={100*(a2-a1)/a1:+.2f}%")

    # front intersections under the same criterion
    def cut(field, nodes_src, tris_src, x0, ys, quadratic):
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
            if quadratic:
                vals.append(float(np.dot(field[tris6[e]], A.shape_vals(l1[j], l2[j]))))
            else:
                vals.append(float(l1[j]*field[tris1[e, 0]] +
                                  l2[j]*field[tris1[e, 1]] +
                                  (1-l1[j]-l2[j])*field[tris1[e, 2]]))
        return np.array(vals)

    ys = np.arange(-0.6, 0.6001, 0.002)
    eps = 1e-8
    print("\n[front cuts] common criterion eps=1e-8")
    cuts = {}
    for x0 in (2.35, 2.38, 2.39, 2.40, 2.41):
        v1 = cut(u1, p1["nodes"], tris1, x0, ys, False)
        v2 = cut(u2, nodes6, tris6, x0, ys, True)
        def intervals(v):
            neu = (~np.isnan(v)) & (np.abs(v) <= eps)
            if not neu.any():
                return []
            idx = np.flatnonzero(neu); res = []
            s = idx[0]; p = idx[0]
            for k in idx[1:]:
                if k != p+1:
                    res.append((float(ys[s]), float(ys[p]))); s = k
                p = k
            res.append((float(ys[s]), float(ys[p])))
            return res
        cuts[x0] = dict(p1=intervals(v1), p2=intervals(v2))
        print(f"  x={x0:.2f}: P1={[(round(a,3), round(b,3)) for a,b in cuts[x0]['p1']]}")
        print(f"           P2={[(round(a,3), round(b,3)) for a,b in cuts[x0]['p2']]}")

    json.dump(dict(areas={str(k): v for k, v in out.items()}, cuts=cuts),
              open("p2_common_zero.json", "w"), indent=2)
    print("\nwrote p2_common_zero.json")


if __name__ == "__main__":
    main()
