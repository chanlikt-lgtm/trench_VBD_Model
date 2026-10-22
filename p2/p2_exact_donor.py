#!/usr/bin/env python3
"""Exact partial-element donor load for the P2 obstacle problem.

f_i = rho * int_{T cap {x >= xb}} N_i(x,y) dA

Implemented entirely in the parent triangle's barycentric frame:
  * clip T against x >= xb  -> convex polygon (<= 4 vertices, physical);
  * fan-triangulate the polygon;
  * for each sub-triangle, map its vertices to parent barycentrics and apply a
    3-point degree-2 exact rule (weights 1/3), with sub-area from the
    barycentric determinant scaled by the parent area.

Sanity identities used as tests:
  * unclipped element: int N = [0,0,0,A/3,A/3,A/3];
  * sum_i int N_i over any clipped piece = its area;
  * sum f = rho * exact charged area.
"""
from __future__ import annotations
import numpy as np

import p2_audit as A
from p2_bernstein import B  # noqa: F401  (kept for import parity)
from p2_charge_mask_audit import clip_area

M = A.M


def clip_poly(P, xb):
    """Clip a convex polygon (list of (x,y)) to x >= xb."""
    out = []
    n = len(P)
    for i in range(n):
        cur = P[i]; nxt = P[(i+1) % n]
        cin = cur[0] >= xb; nin = nxt[0] >= xb
        if cin:
            out.append(cur)
        if cin != nin:
            t = (xb - cur[0])/(nxt[0] - cur[0])
            out.append((xb, cur[1] + t*(nxt[1] - cur[1])))
    return out


def bary(P, pt):
    Mx = np.array([[P[1][0]-P[0][0], P[2][0]-P[0][0]],
                   [P[1][1]-P[0][1], P[2][1]-P[0][1]]])
    a, b = np.linalg.solve(Mx, np.array(pt) - np.array(P[0]))
    return (1.0 - a - b, a, b)      # weights of v0, v1, v2


def int_N(P, xb, area):
    """Integrals of the six P2 shape functions over T cap {x >= xb}."""
    poly = clip_poly([tuple(p) for p in P], xb)
    if len(poly) < 3:
        return np.zeros(6)
    bw = [bary(P, p) for p in poly]
    acc = np.zeros(6)
    W = [(2/3, 1/6, 1/6), (1/6, 2/3, 1/6), (1/6, 1/6, 2/3)]
    for k in range(1, len(bw)-1):
        tri = (bw[0], bw[k], bw[k+1])
        d = abs((tri[1][0]-tri[0][0])*(tri[2][1]-tri[0][1]) -
                (tri[2][0]-tri[0][0])*(tri[1][1]-tri[0][1]))
        sub_area = d*area
        for w in W:
            L1 = sum(wj*bj[0] for wj, bj in zip(w, tri))
            L2 = sum(wj*bj[1] for wj, bj in zip(w, tri))
            L3 = sum(wj*bj[2] for wj, bj in zip(w, tri))
            acc += (sub_area/3.0) * np.array(
                [L1*(2*L1-1), L2*(2*L2-1), L3*(2*L3-1),
                 4*L1*L2, 4*L2*L3, 4*L3*L1])
    return acc


def exact_load(nodes, tris, mat, xb=None):
    xb = A.BODY_BOTTOM if xb is None else xb
    area = np.abs(A.element_geom(nodes, tris)[2])
    P = nodes[tris][:, :3]
    xmax = P[:, :, 0].max(1)
    si = mat == M.SI
    cand = np.flatnonzero(si & (xmax > xb))
    f = np.zeros(len(nodes))
    total = 0.0
    for e in cand:
        I = int_N(P[e], xb, area[e])
        f[tris[e]] += A.RHO * I
        total += I.sum()
    return f, float(total)


def test(nodes, tris, mat, tag):
    area = np.abs(A.element_geom(nodes, tris)[2])
    P = nodes[tris][:, :3]
    si = mat == M.SI
    full = np.flatnonzero(si & (P[:, :, 0].min(1) > A.BODY_BOTTOM))
    e = int(full[0])
    I = int_N(P[e], -1e9, area[e])
    err = np.abs(I - np.array([0, 0, 0, area[e]/3]*3 if False else
                              [0, 0, 0, area[e]/3, area[e]/3, area[e]/3])).max()
    f, tot = exact_load(nodes, tris, mat)
    print(f"[{tag}] full-element identity max err = {err:.2e}")
    print(f"[{tag}] sum f / rho = {tot:.12f}")
    return err, tot, f


def main():
    for tag, fn in (("B0", "p2_B0_p2mesh.npz"), ("B1", "p2_B1_p2mesh.npz")):
        z = np.load(fn)
        test(z["nodes"], z["tris"], z["mat"], tag)


if __name__ == "__main__":
    main()
