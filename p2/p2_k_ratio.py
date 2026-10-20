#!/usr/bin/env python3
"""Is T'K2T a uniform multiple of K1? If the P2 quadrature weight is doubled,
the ratio is exactly 2 and the mismatch norm equals ||K1||.
"""
from __future__ import annotations
import numpy as np
import scipy.sparse as sp

import p2_audit as A
from p2_nested_consistency import build_T, assemble_p1, assemble_p2

M = A.M


def main():
    z = np.load("p2_r3_solution.npz")
    nodes6, tris6, mat = z["nodes"], z["tris"], z["mat"]
    T, o2n, verts = build_T(nodes6, tris6)
    nodes = nodes6[verts]
    tris = o2n[tris6[:, :3]]
    cxs = nodes6[tris6][:, :3, 0].mean(axis=1)
    si = mat == M.SI
    body = si & (cxs <= A.BODY_BOTTOM + 1e-12)
    charge = si & ~body

    K2, _ = assemble_p2(nodes6, tris6, mat, charge)
    K1, _ = assemble_p1(nodes, tris, mat, charge)
    K2T = T.T @ (K2 @ T)

    for c in (1.0, 2.0, 0.5):
        D = K2T - c*K1
        print(f"c={c}: ||K2T - c K1||_F / ||c K1||_F = "
              f"{np.sqrt((D.multiply(D)).sum())/np.sqrt((K2T.multiply(K2T)).sum()):.6e}")

    d1 = K1.diagonal(); d2 = K2T.diagonal()
    m = (np.abs(d1) > 1e-12)
    r = d2[m]/d1[m]
    print(f"diag ratio: median={np.median(r):.9f}  p05={np.percentile(r,5):.9f} "
          f"p95={np.percentile(r,95):.9f}")

    # direct quadrature check: integral of 1 over one element
    print("\n[quadrature] sum(QW) =", A.QW.sum(), " |det| = 2A")
    print("unit integral via rule = sum(QW)*|det| =", A.QW.sum()*2.0, "* A")
    print("reference-triangle area = 0.5, so weights should sum to 0.5")


if __name__ == "__main__":
    main()
