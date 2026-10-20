#!/usr/bin/env python3
"""P2.1 diagnostics — measure, do not infer.

A. print the six raw shape integrals / A and the six load entries for 3 Si elements
B. rebuild adjacency from ALL SIX DOFs and recount anchored components
C. inspect K_FF rows directly (uncoupled rows, symmetry)
D. corrected manufactured-solution tests (zero-source and nonzero-source)
"""
from __future__ import annotations
import json
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.sparse.csgraph import connected_components

import p2_audit as A

M = A.M
QP, QW = A.QP, A.QW
RHO, EPS_SI, EPS_OX = A.RHO, A.EPS_SI, A.EPS_OX
V_DRAIN = A.V_DRAIN
out = {}


# ---------- A. raw integrals ----------
def diag_a(nodes, tris, mat):
    p, det, area, g1, g2, g3 = A.element_geom(nodes, tris)
    si = np.flatnonzero(mat == M.SI)
    rows = []
    for m in si[[0, len(si)//2, len(si)-1]]:
        I = np.zeros(6)
        for q in range(3):
            L1, L2 = QP[q]
            I += QW[q] * abs(det[m]) * A.shape_vals(L1, L2)
        Ar = abs(area[m])
        rows.append(dict(element=int(m), area=float(Ar),
                         over_area=[float(v) for v in I/Ar],
                         load_over_rhoA=[float(v) for v in I/Ar]))
    return rows


# ---------- B. full six-DOF adjacency ----------
def diag_b(nodes, tris):
    n = len(nodes)
    r, c = [], []
    for e in tris:
        for i in range(6):
            for j in range(6):
                if i != j:
                    r.append(e[i]); c.append(e[j])
    Adj = sp.coo_matrix((np.ones(len(r)), (r, c)), shape=(n, n)).tocsr()
    return connected_components(Adj, directed=False)


def main():
    d = np.load("p2_r3_p2mesh.npz")
    nodes, tris, mat = d["nodes"], d["tris"], d["mat"]

    out["A_raw_integrals"] = diag_a(nodes, tris, mat)
    print("[A] six raw integrals / A  (expect 0 0 0 1/3 1/3 1/3)")
    for r in out["A_raw_integrals"]:
        s = " ".join(f"{v:+.6f}" for v in r["over_area"])
        print(f"    el {r['element']:>6}  A={r['area']:.3e}  [{s}]")

    # compaction
    used = np.unique(tris.ravel())
    o2n = np.full(len(nodes), -1, dtype=int)
    o2n[used] = np.arange(len(used))
    nodes_c = nodes[used]; tris_c = o2n[tris]
    n = len(nodes_c)

    K, f = A.assemble(nodes_c, tris_c, mat)
    x, y = nodes_c[:, 0], nodes_c[:, 1]
    nodes_of = {}
    for code in (M.SI, M.OX, M.GATE, M.SHIELD):
        sel = mat == code
        nodes_of[code] = set(np.unique(tris_c[sel]).tolist()) if sel.any() else set()
    si_n, ox_n = nodes_of[M.SI], nodes_of[M.OX]
    gate_n, shield_n = nodes_of[M.GATE], nodes_of[M.SHIELD]
    src = {i for i in si_n if x[i] <= 1e-12}
    drn = {i for i in si_n if x[i] >= M.BOTTOM_DEPTH - 1e-9}
    body = {i for i in si_n if x[i] <= A.BODY_BOTTOM + 1e-12}

    fixed = np.zeros(n, bool); val = np.zeros(n)
    for s, v in ((src, A.V_SRC), (body, A.V_SRC), (drn, V_DRAIN),
                 (gate_n, A.V_G), (shield_n, A.V_SG)):
        for i in s:
            fixed[i] = True; val[i] = v
    idx = np.flatnonzero(~fixed)

    ncomp, lab = diag_b(nodes_c, tris_c)
    unanch = [c for c in range(ncomp) if not ((lab == c) & fixed).any()]
    sizes = sorted([int((lab == c).sum()) for c in range(ncomp)], reverse=True)
    out["B_components"] = dict(components=int(ncomp), unanchored=int(len(unanch)),
                               largest=sizes[:5], anchored=int(ncomp-len(unanch)))
    print(f"[B] full six-DOF graph: components={ncomp} unanchored={len(unanch)} "
          f"largest={sizes[:5]}")

    dg = K.diagonal()
    offdiag = np.array([K[i].nnz - 1 for i in idx])
    sym = abs(K - K.T).max() if K.shape[0] else 0.0
    out["C_matrix"] = dict(zero_free_diag=int((dg[idx] <= 0).sum()),
                           min_free_diag=float(dg[idx].min()),
                           rows_with_no_offdiag=int((offdiag == 0).sum()),
                           max_asymmetry=float(sym),
                           nnz=int(K.nnz))
    print(f"[C] K_FF: zeroDiag={(dg[idx]<=0).sum()} rowsNoOffdiag={(offdiag==0).sum()} "
          f"maxAsym={sym:.2e}")

    # ---------- D. manufactured tests ----------
    # D1: phi = x^2 - y^2, zero source, constant eps, Dirichlet everywhere outside
    eps_t = 1.0
    def assemble_const(nodes, tris, eps):
        p, det, area, g1, g2, g3 = A.element_geom(nodes, tris)
        r = np.zeros((len(tris), 36), dtype=np.int64)
        c = np.zeros((len(tris), 36), dtype=np.int64)
        v = np.zeros((len(tris), 36))
        for m in range(len(tris)):
            ke = np.zeros((6, 6))
            for q in range(3):
                L1, L2 = QP[q]
                G = A.shape_grads(L1, L2, g1[m], g2[m], g3[m])
                ke += QW[q] * abs(det[m]) * eps * (G @ G.T)
            ids = tris[m]
            r[m] = np.repeat(ids, 6); c[m] = np.tile(ids, 6); v[m] = ke.ravel()
        return sp.coo_matrix((v.ravel(), (r.ravel(), c.ravel())), shape=(len(nodes),)*2).tocsr()

    # use a small patch of the real mesh but treat it as single-material/constant-eps
    Kc = assemble_const(nodes_c, tris_c, eps_t)
    ph1 = x**2 - y**2
    # anchor: everything on the outer box is Dirichlet with the exact value
    outer = (x <= 1e-12) | (x >= M.BOTTOM_DEPTH-1e-9) | (np.abs(y) >= M.Y_EDGE-1e-9)
    free1 = ~outer
    res1 = np.zeros(n)
    res1[free1] = (Kc @ ph1)[free1]
    scale1 = np.abs(Kc @ ph1).max() + 1e-30
    out["D1_zero_source"] = dict(max_interior=float(np.abs(res1[free1]).max()),
                                 relative=float(np.abs(res1[free1]).max()/scale1))
    print(f"[D1] phi=x^2-y^2 zero source: interior residual "
          f"{np.abs(res1[free1]).max():.3e} (rel {np.abs(res1[free1]).max()/scale1:.2e})")

    # D2: phi = x^2 + y^2, source -4 eps, consistent P2 load on midside nodes
    ph2 = x**2 + y**2
    load = np.zeros(n)
    p, det, area, g1, g2, g3 = A.element_geom(nodes_c, tris_c)
    for m in range(len(tris_c)):
        I = np.zeros(6)
        for q in range(3):
            L1, L2 = QP[q]
            I += QW[q] * abs(det[m]) * A.shape_vals(L1, L2)
        load[tris_c[m]] += 4.0 * eps_t * I
    res2 = np.zeros(n)
    res2[free1] = (Kc @ ph2 + load)[free1]
    scale2 = np.abs(Kc @ ph2).max() + np.abs(load).max() + 1e-30
    out["D2_nonzero_source"] = dict(max_interior=float(np.abs(res2[free1]).max()),
                                    relative=float(np.abs(res2[free1]).max()/scale2))
    print(f"[D2] phi=x^2+y^2 source -4eps: interior residual "
          f"{np.abs(res2[free1]).max():.3e} (rel {np.abs(res2[free1]).max()/scale2:.2e})")

    with open("p2_diag.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print("wrote p2_diag.json")


if __name__ == "__main__":
    main()
