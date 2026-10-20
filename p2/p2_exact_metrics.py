#!/usr/bin/env python3
"""Controlled-error neutral area and analytic front cuts (corrected).

P1: exact triangle clipping by u <= eps.
P2: equal-area sub-triangle quadrature, K = 8, 16, 32.
Fronts: exact element/line intersection, u(t) sampled on the segment with
bisection at every crossing; intervals merged with a strict geometric union.
"""
from __future__ import annotations
import json
import numpy as np

import p2_audit as A

M = A.M


def sub_points(K):
    L1, L2 = [], []
    for i in range(K):
        for j in range(K-i):
            L1.append((i+1/3)/K); L2.append((j+1/3)/K)
    return np.array(L1), np.array(L2)


def area_p2(nodes, tris, mat, u, eps, K):
    area = np.abs(A.element_geom(nodes, tris)[2])
    cxs = nodes[tris][:, :3, 0].mean(axis=1)
    charged = (mat == M.SI) & (cxs > A.BODY_BOTTOM + 1e-12)
    els = np.flatnonzero(charged)
    L1, L2 = sub_points(K)
    N = np.array([A.shape_vals(a, b) for a, b in zip(L1, L2)])
    tot = 0.0
    for s in range(0, len(els), 1000):
        b = els[s:s+1000]
        U = u[tris[b]] @ N.T
        tot += float(((U <= eps).mean(axis=1)*area[b]).sum())
    return tot


def area_p1(nodes, tris, mat, u, eps):
    area = np.abs(A.element_geom(nodes, tris)[2])
    cxs = nodes[tris[:, :3], 0].mean(axis=1)
    charged = (mat == M.SI) & (cxs > A.BODY_BOTTOM + 1e-12)
    els = np.flatnonzero(charged)
    tot = 0.0
    for e in els:
        ids = tris[e, :3]
        v = u[ids]
        if np.all(v <= eps):
            tot += area[e]; continue
        if np.all(v > eps):
            continue
        poly = []
        P = nodes[ids]
        for k in range(3):
            p, q = P[k], P[(k+1) % 3]
            vp, vq = v[k], v[(k+1) % 3]
            if vp <= eps:
                poly.append(p)
            if (vp-eps)*(vq-eps) < 0:
                t = (eps-vp)/(vq-vp)
                poly.append(p + t*(q-p))
        if len(poly) >= 3:
            x = np.array([p[0] for p in poly]); y = np.array([p[1] for p in poly])
            tot += 0.5*abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
    return tot


def cut_intervals(nodes, tris, mat, u, x0, eps, quad):
    ids3 = tris[:, :3]
    P = nodes[ids3]
    xmin = P[:, :, 0].min(axis=1); xmax = P[:, :, 0].max(axis=1)
    cxs = P[:, :, 0].mean(axis=1)
    charged = (mat == M.SI) & (cxs > A.BODY_BOTTOM + 1e-12)
    cand = np.flatnonzero(charged & (xmin-1e-12 <= x0) & (xmax+1e-12 >= x0))
    Ns = 65
    segs = []
    for e in cand:
        pts = []
        for k in range(3):
            p, q = P[e, k], P[e, (k+1) % 3]
            if abs(p[0]-q[0]) < 1e-15:
                if abs(p[0]-x0) < 1e-12:
                    pts.append(p.copy()); pts.append(q.copy())
                continue
            t = (x0-p[0])/(q[0]-p[0])
            if 0.0 <= t <= 1.0:
                pts.append(p + t*(q-p))
        if len(pts) < 2:
            continue
        yy = sorted({round(float(p[1]), 14) for p in pts})
        if len(yy) < 2:
            continue
        y0, y1 = yy[0], yy[-1]
        if abs(y1-y0) < 1e-12:
            continue
        Mx = np.array([[P[e, 1, 0]-P[e, 0, 0], P[e, 2, 0]-P[e, 0, 0]],
                       [P[e, 1, 1]-P[e, 0, 1], P[e, 2, 1]-P[e, 0, 1]]])
        Mxinv = np.linalg.inv(Mx)

        def val_at(t):
            st = Mxinv @ (np.array([x0, y0 + t*(y1-y0)]) - P[e, 0])
            l1, l2 = st[0], st[1]
            if quad:
                return float(np.dot(u[tris[e]], A.shape_vals(l1, l2)))
            return float(l1*u[tris[e, 0]] + l2*u[tris[e, 1]] +
                         (1-l1-l2)*u[tris[e, 2]])

        ts = np.linspace(0, 1, Ns)
        vals = np.array([val_at(t) for t in ts])
        ind = vals <= eps
        i = 0
        while i < Ns:
            if not ind[i]:
                i += 1; continue
            j = i
            while j+1 < Ns and ind[j+1]:
                j += 1
            if i == 0:
                t_start = 0.0
            else:
                lo, hi = ts[i-1], ts[i]
                for _ in range(60):
                    mid = 0.5*(lo+hi)
                    if val_at(mid) <= eps: hi = mid
                    else: lo = mid
                t_start = 0.5*(lo+hi)
            if j == Ns-1:
                t_end = 1.0
            else:
                lo, hi = ts[j], ts[j+1]
                for _ in range(60):
                    mid = 0.5*(lo+hi)
                    if val_at(mid) <= eps: lo = mid
                    else: hi = mid
                t_end = 0.5*(lo+hi)
            segs.append((y0 + t_start*(y1-y0), y0 + t_end*(y1-y0)))
            i = j+1
    segs.sort()
    merged = []
    for a, b in segs:
        if b-a < 1e-9:
            continue
        if merged and a <= merged[-1][1] + 2e-4:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    return merged


def main():
    eps = 1e-8
    out = {}
    print(f"eps = {eps:g}")
    p2 = {}
    for tag, fn in (("r3", "p2_p3_r3_polished.npz"), ("r4", "p2_p3_r4_polished.npz"),
                    ("mid", "p2_p3_rmid_polished.npz")):
        z = np.load(fn)
        nodes, tris, mat, u = z["nodes"], z["tris"], z["mat"], z["u"]
        a8 = area_p2(nodes, tris, mat, u, eps, 8)
        a16 = area_p2(nodes, tris, mat, u, eps, 16)
        a32 = area_p2(nodes, tris, mat, u, eps, 32)
        cuts = {x: cut_intervals(nodes, tris, mat, u, x, eps, True)
                for x in (2.35, 2.38, 2.39, 2.40, 2.41)}
        p2[tag] = dict(a8=a8, a16=a16, a32=a32, cuts=cuts)
        print(f"[P2 {tag}] K=8 {a8:.9f}  K=16 {a16:.9f}  K=32 {a32:.9f}  "
              f"(K16-K8 {a16-a8:+.3e}, K32-K16 {a32-a16:+.3e})")
        for x, iv in cuts.items():
            print(f"   x={x:.2f}: {[(round(a, 4), round(b, 4)) for a, b in iv]}")

    p1 = {}
    for tag, fn in (("r3", "p1_on_p2mesh_obstacle.npz"),
                    ("r4", "p1_on_p2mesh_obstacle_r4.npz")):
        z = np.load(fn)
        nodes, tris, mat, u = z["nodes"], z["tris"], z["mat"], z["u"]
        a = area_p1(nodes, tris, mat, u, eps)
        cuts = {x: cut_intervals(nodes, tris, mat, u, x, eps, False)
                for x in (2.35, 2.38, 2.39, 2.40, 2.41)}
        p1[tag] = dict(a=a, cuts=cuts)
        print(f"[P1 {tag}] exact area {a:.9f}")
        for x, iv in cuts.items():
            print(f"   x={x:.2f}: {[(round(a_, 4), round(b_, 4)) for a_, b_ in iv]}")

    print("\n[ratios]")
    for tag in ("r3", "r4"):
        print(f"  {tag}: P1 exact {p1[tag]['a']:.9f}  P2 K32 {p2[tag]['a32']:.9f}  "
              f"P2/P1 {p2[tag]['a32']/p1[tag]['a']:.5f}")

    json.dump(dict(p2=p2, p1=p1), open("p2_exact_metrics.json", "w"), indent=2)
    print("\nwrote p2_exact_metrics.json")


if __name__ == "__main__":
    main()
