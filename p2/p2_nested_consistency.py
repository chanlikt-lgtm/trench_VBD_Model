#!/usr/bin/env python3
"""Nested-space consistency: P1 embedded in P2.

Prolongation T (n_P2 x n_P1):
  vertex  -> identity
  midside -> 0.5 * (two endpoints)

Checks, before Dirichlet elimination, identical geometry/materials/charge:
    K_P1  vs  T' K_P2 T
    f_P1  vs  T' f_P2
    T' T  vs  I

Run on any P2 mesh that has a solution archive (r3 or r4).
"""
from __future__ import annotations
import sys
import numpy as np
import scipy.sparse as sp

import p2_audit as A

M = A.M
QP, QW = A.QP, A.QW


def geom(nodes, tris):
    p = nodes[tris]
    x1, y1 = p[:, 0, 0], p[:, 0, 1]
    x2, y2 = p[:, 1, 0], p[:, 1, 1]
    x3, y3 = p[:, 2, 0], p[:, 2, 1]
    det = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)
    g1 = np.stack([(y2-y3), (x3-x2)], 1)/det[:, None]
    g2 = np.stack([(y3-y1), (x1-x3)], 1)/det[:, None]
    g3 = np.stack([(y1-y2), (x2-x1)], 1)/det[:, None]
    return det, 0.5*np.abs(det), g1, g2, g3


def assemble_p2(nodes, tris, mat, charge):
    n = len(nodes)
    det, area, g1, g2, g3 = geom(nodes, tris)
    eps = np.where(mat == M.OX, A.EPS_OX, A.EPS_SI)
    G = np.zeros((len(tris), 3, 6, 2))
    for q in range(3):
        L1, L2 = QP[q]
        L3 = 1.0 - L1 - L2
        G[:, q, 0] = (4*L1-1)*g1
        G[:, q, 1] = (4*L2-1)*g2
        G[:, q, 2] = (4*L3-1)*g3
        G[:, q, 3] = 4*(L1*g2 + L2*g1)
        G[:, q, 4] = 4*(L2*g3 + L3*g2)
        G[:, q, 5] = 4*(L3*g1 + L1*g3)
    ke = np.zeros((len(tris), 6, 6))
    for q in range(3):
        w = QW[q]*np.abs(det)*eps
        ke += w[:, None, None]*np.einsum('eip,ejp->eij', G[:, q], G[:, q])
    ids = tris
    rows = np.repeat(ids, 6, axis=1).ravel()
    cols = np.tile(ids, (1, 6)).ravel()
    K = sp.coo_matrix((ke.reshape(-1), (rows, cols)), shape=(n, n)).tocsr()
    f = np.zeros(n)
    for e in np.flatnonzero(charge):
        f[tris[e][3:]] += A.RHO*area[e]/3.0
    return K, f


def assemble_p1(nodes, tris, mat, charge):
    n = len(nodes)
    p = nodes[tris]
    det = ((p[:, 1, 0]-p[:, 0, 0])*(p[:, 2, 1]-p[:, 0, 1]) -
           (p[:, 2, 0]-p[:, 0, 0])*(p[:, 1, 1]-p[:, 0, 1]))
    area = 0.5*np.abs(det)
    b = np.stack([p[:, 1, 1]-p[:, 2, 1], p[:, 2, 1]-p[:, 0, 1], p[:, 0, 1]-p[:, 1, 1]], 1)
    c = np.stack([p[:, 2, 0]-p[:, 1, 0], p[:, 0, 0]-p[:, 2, 0], p[:, 1, 0]-p[:, 0, 0]], 1)
    eps = np.where(mat == M.OX, A.EPS_OX, A.EPS_SI)
    ke = (eps/(4*area))[:, None, None]*(
        np.einsum('ei,ej->eij', b, b) + np.einsum('ei,ej->eij', c, c))
    ids = tris
    rows = np.repeat(ids, 3, axis=1).ravel()
    cols = np.tile(ids, (1, 3)).ravel()
    K = sp.coo_matrix((ke.reshape(-1), (rows, cols)), shape=(n, n)).tocsr()
    f = np.zeros(n)
    for e in np.flatnonzero(charge):
        f[tris[e]] += A.RHO*area[e]/3.0
    return K, f


def build_T(nodes6, tris6):
    n6 = len(nodes6)
    verts = np.unique(tris6[:, :3])
    o2n = np.full(n6, -1, dtype=int)
    o2n[verts] = np.arange(len(verts))
    rows = list(verts)
    cols = list(o2n[verts])
    vals = [1.0]*len(verts)
    seen = {}
    for t in tris6:
        for i, j, k in ((t[0], t[1], t[3]), (t[1], t[2], t[4]), (t[2], t[0], t[5])):
            if k not in seen:
                seen[k] = (i, j)
    for k, (i, j) in seen.items():
        rows += [k, k]
        cols += [o2n[i], o2n[j]]
        vals += [0.5, 0.5]
    return sp.coo_matrix((vals, (rows, cols)), shape=(n6, len(verts))).tocsr(), o2n, verts


def run(fn):
    z = np.load(fn)
    nodes6, tris6, mat = z["nodes"], z["tris"], z["mat"]
    T, o2n, verts = build_T(nodes6, tris6)
    nodes = nodes6[verts]
    tris = o2n[tris6[:, :3]]

    cxs = nodes6[tris6][:, :3, 0].mean(axis=1)
    si = mat == M.SI
    body = si & (cxs <= A.BODY_BOTTOM + 1e-12)
    charge = si & ~body

    K2, f2 = assemble_p2(nodes6, tris6, mat, charge)
    K1, f1 = assemble_p1(nodes, tris, mat, charge)

    TtK2T = T.T @ (K2 @ T)
    D = K1 - TtK2T
    nK = np.sqrt((K1.multiply(K1)).sum())
    print(f"\n{fn}: n1={len(nodes)} n2={len(nodes6)} elems={len(tris6)}")
    print(f"  ||K1 - T'K2T||_F / ||K1||_F = {np.sqrt((D.multiply(D)).sum())/nK:.3e}")
    print(f"  max |K1 - T'K2T|             = {np.abs(D).max():.3e}")
    df = f1 - T.T @ f2
    print(f"  ||f1 - T'f2|| / ||f1||       = {np.linalg.norm(df)/np.linalg.norm(f1):.3e}")
    print(f"  max |f1 - T'f2|              = {np.abs(df).max():.3e}")
    TT = T.T @ T
    I = sp.eye(TT.shape[0], format="csr")
    E = TT - I
    print(f"  ||T'T - I||_F                = {np.sqrt((E.multiply(E)).sum()):.3e}")
    return dict(K_rel=float(np.sqrt((D.multiply(D)).sum())/nK),
                K_max=float(np.abs(D).max()),
                f_rel=float(np.linalg.norm(df)/np.linalg.norm(f1)),
                f_max=float(np.abs(df).max()))


if __name__ == "__main__":
    import json
    out = {}
    for tag, fn in (("r3", "p2_r3_solution.npz"), ("r4", "p2_r4_solution.npz")):
        out[tag] = run(fn)
    json.dump(out, open("p2_nested_consistency.json", "w"), indent=2)
    print("\nwrote p2_nested_consistency.json")
