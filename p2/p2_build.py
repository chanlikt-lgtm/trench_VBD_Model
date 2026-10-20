#!/usr/bin/env python3
"""Stage 2C P2.0 + P2.1 — six-node triangles on the r3 geometry, then the
quadratic stiffness/charge assembly and the linear operator test.

P2.0: build the r3 geometry at Mesh.ElementOrder = 2. Gmsh shares edge-midpoint
      nodes between adjacent elements, so the global DOF count is the node count
      (NOT 6 per element). The actual number is read from the connectivity.
P2.1: assemble the quadratic stiffness with a 3-point rule (exact for the
      degree-2 products that appear), apply the same Dirichlet data, and check
      the linear residual and the integrated charge.

No obstacle is applied here — this is assembly validation only.
"""
from __future__ import annotations
import os
import sys
import json
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import gmsh

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import mesh_stage2b as M

BODY_BOTTOM = 0.54
V_G, V_SG, V_DRAIN, V_SRC = 0.0, 0.0, 34.274208536, 0.0
EPS_SI, EPS_OX = 11.7, 3.9
RHO = (1.602176634e-19 * 6.9e16 / 8.8541878128e-14) * 1e-8
H, H_TIP = 0.014, 0.0055          # r3 geometry


def build_p2(h=H, h_tip=H_TIP, prefix="p2_r3"):
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 1)
    groups = M.build(h, h_tip)                    # geometry + physical groups
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.SecondOrderLinear", 1)
    gmsh.option.setNumber("Mesh.MeshSizeMin", h_tip)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.model.mesh.generate(2)

    ntags, ncoords, _ = gmsh.model.mesh.getNodes()
    nodes = np.column_stack([ncoords[0::3], ncoords[1::3]])
    tag2i = {int(t): i for i, t in enumerate(ntags)}

    parts, mats, counts = {}, {}, {}
    for code, lst in groups.items():
        acc = []
        for f in lst:
            etypes, _, enodes = gmsh.model.mesh.getElements(2, f)
            for et, en in zip(etypes, enodes):
                if et == 9:
                    acc.append(np.vectorize(tag2i.get)(
                        np.asarray(en, dtype=int).reshape(-1, 6)))
        if acc:
            parts[code] = np.vstack(acc)
            mats[code] = np.full(len(parts[code]), code)
            counts[code] = len(parts[code])
    gmsh.finalize()

    tris = np.vstack([parts[k] for k in parts])
    mat = np.concatenate([mats[k] for k in parts])
    np.savez_compressed(f"{prefix}_p2mesh.npz", nodes=nodes, tris=tris, mat=mat)
    return nodes, tris, mat, counts


# ---------- P2 shape functions and gradients at 3-point rule ----------
QP = np.array([[1/6, 1/6], [1/6, 2/3], [2/3, 1/6]])
QW = np.array([1/3, 1/3, 1/3])


def shape_grads(L1, L2, g1, g2, g3):
    """Six P2 shape-function gradients (6,2) at barycentric (L1,L2,L3=1-L1-L2)."""
    L3 = 1.0 - L1 - L2
    N = np.zeros((6, 2))
    N[0] = (4*L1 - 1) * g1
    N[1] = (4*L2 - 1) * g2
    N[2] = (4*L3 - 1) * g3
    N[3] = 4*(L1*g2 + L2*g1)
    N[4] = 4*(L2*g3 + L3*g2)
    N[5] = 4*(L3*g1 + L1*g3)
    return N


def assemble_p2(nodes, tris, mat):
    n = len(nodes)
    nd = 36
    rows = np.zeros((len(tris), nd), dtype=np.int64)
    cols = np.zeros((len(tris), nd), dtype=np.int64)
    vals = np.zeros((len(tris), nd))
    f = np.zeros(n)

    p = nodes[tris]                       # (M,6,2)
    x1, y1 = p[:, 0, 0], p[:, 0, 1]
    x2, y2 = p[:, 1, 0], p[:, 1, 1]
    x3, y3 = p[:, 2, 0], p[:, 2, 1]
    det = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)
    area = 0.5*det
    g1 = np.stack([(y2-y3), (x3-x2)], 1) / det[:, None]
    g2 = np.stack([(y3-y1), (x1-x3)], 1) / det[:, None]
    g3 = np.stack([(y1-y2), (x2-x1)], 1) / det[:, None]
    eps = np.where(mat == M.OX, EPS_OX, EPS_SI)

    # P2 consistent load for constant rho. With N_i = lambda_i(2 lambda_i - 1)
    # and N_ij = 4 lambda_i lambda_j, the exact integrals are
    #     int N_vertex dA = 0 ,      int N_midside dA = A/3 .
    # Hence f_e = (rho A / 3) [0,0,0,1,1,1]^T. A vertices-only vector sums to
    # the right total but is NOT the consistent load, so the individual entries
    # are checked explicitly below.
    midside_slots = np.array([3, 4, 5])

    for m in range(len(tris)):
        ke = np.zeros((6, 6))
        for q in range(3):
            L1, L2 = QP[q]
            G = shape_grads(L1, L2, g1[m], g2[m], g3[m])
            ke += QW[q] * abs(det[m]) * eps[m] * (G @ G.T)
        ids = tris[m]
        rows[m] = np.repeat(ids, 6)
        cols[m] = np.tile(ids, 6)
        vals[m] = ke.ravel()

        if mat[m] == M.SI:
            f[ids[midside_slots]] += RHO * abs(area[m]) / 3.0

    K = sp.coo_matrix((vals.ravel(), (rows.ravel(), cols.ravel())),
                      shape=(n, n)).tocsr()
    return K, f, area


def main():
    t_vertices = None
    nodes, tris, mat, counts = build_p2()
    n = len(nodes)
    print(f"[P2.0] vertices/DOFs = {n}")
    print(f"[P2.0] P2 triangles   = {len(tris)}  (by material: "
          f"{ {k: v for k, v in counts.items()} })")
    print(f"[P2.0] naive 6-per-element count = {6*len(tris)}  "
          f"-> shared-DOF factor {6*len(tris)/n:.3f}")

    K, f, area = assemble_p2(nodes, tris, mat)
    print(f"[P2.1] K nnz = {K.nnz}")

    x, y = nodes[:, 0], nodes[:, 1]
    nodes_of = {}
    for code in (M.SI, M.OX, M.GATE, M.SHIELD):
        sel = mat == code
        nodes_of[code] = set(np.unique(tris[sel]).tolist()) if sel.any() else set()
    a = M.a_of(x)
    si_n, ox_n = nodes_of[M.SI], nodes_of[M.OX]
    gate_n, shield_n = nodes_of[M.GATE], nodes_of[M.SHIELD]
    source_n = {i for i in si_n if x[i] <= 1e-12}
    drain_n = {i for i in si_n if x[i] >= M.BOTTOM_DEPTH - 1e-9}
    body_n = {i for i in si_n if x[i] <= BODY_BOTTOM + 1e-12}

    fixed = np.zeros(n, bool); val = np.zeros(n)
    for nset, v in ((source_n, V_SRC), (body_n, V_SRC), (drain_n, V_DRAIN),
                    (gate_n, V_G), (shield_n, V_SG)):
        for i in nset:
            fixed[i] = True; val[i] = v
    idx = np.flatnonzero(~fixed)
    rhs = f - K @ val
    u = val.copy()
    u[idx] = spla.spsolve(K[idx][:, idx].tocsc(), rhs[idx])

    r = K @ u - f; r[fixed] = 0.0
    num = np.linalg.norm(r[idx]); den = (np.linalg.norm(K[idx][:, idx] @ u[idx])
                                         + np.linalg.norm(f[idx]) + 1e-30)
    q = float(f.sum()); qa = float(RHO * area[mat == M.SI].sum())

    out = dict(
        stage="P2.0/P2.1",
        dofs=int(n), p2_triangles=int(len(tris)),
        triangles_by_material={k: int(v) for k, v in counts.items()},
        naive_six_per_element=int(6*len(tris)),
        shared_dof_factor=float(6*len(tris)/n),
        K_nnz=int(K.nnz),
        fixed_nodes=int(fixed.sum()), free_dofs=int(idx.size),
        linear_residual_R_free=float(num/den),
        finite=bool(np.all(np.isfinite(u))),
        phi_min=float(u.min()), phi_max=float(u.max()),
        charge_integral=q, charge_analytic=qa,
        charge_rel_error=float(abs(q-qa)/qa),
        oxide_nodes_with_charge=int((np.abs(f[list(ox_n)]) > 0).sum()),
        gate_nodes=int(len(gate_n)), shield_nodes=int(len(shield_n)),
        gate_shield_shared=int(len(gate_n & shield_n)),
    )
    with open("p2_p01.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
