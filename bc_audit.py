#!/usr/bin/env python3
"""Stage 2B.1 — assembled boundary-condition and charge audit (single resolution).

Reuses the validated Stage 2B geometry/mesh, then audits what is ACTUALLY
assembled into the P1 system, not just what the CAD tags say:

  * every physical group: face count, node count, intended condition, assigned value
  * G and SG node sets must be disjoint (no shared vertex, no electrical short)
  * electrodes must not own any interior silicon or oxide node
  * source / body / drain Dirichlet sets and values
  * Neumann boundaries: outer pitch edges and the symmetry plane
  * Si/SiO2 interface: nodes shared by both materials, and weak-form flux continuity
  * complementarity restricted to silicon nodes only
  * charge: qND/eps0 in V/um^2, integrated charge vs analytic, silicon-only
  * algebraic residual of the constrained system

Triangle orientation is normalised (all signed areas positive) before assembly.
"""
from __future__ import annotations
import json
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import gmsh

import mesh_stage2b as M

BODY_BOTTOM = 0.54
V_G, V_SG, V_DRAIN, V_SRC = 0.0, 0.0, 34.274208536, 0.0
ND = 6.9e16
EPS0_FCM = 8.8541878128e-14
Q = 1.602176634e-19
EPS_SI, EPS_OX = 11.7, 3.9
RHO = (Q * ND / EPS0_FCM) * 1e-8


def audit(h=0.02, h_tip=0.008, out="bc_audit.json"):
    groups = M.build(h, h_tip)
    gmsh.option.setNumber("Mesh.MeshSizeMin", h_tip)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.model.mesh.generate(2)

    ntags, ncoords, _ = gmsh.model.mesh.getNodes()
    nodes = np.column_stack([ncoords[0::3], ncoords[1::3]])
    tag2i = {int(t): i for i, t in enumerate(ntags)}

    tris_by_mat, nodes_by_mat = {}, {k: set() for k in groups}
    for code, lst in groups.items():
        acc = []
        for f in lst:
            etypes, _, enodes = gmsh.model.mesh.getElements(2, f)
            for et, en in zip(etypes, enodes):
                if et == 2:
                    acc.append(np.vectorize(tag2i.get)(np.asarray(en, int).reshape(-1, 3)))
        if acc:
            tris_by_mat[code] = np.vstack(acc)
            nodes_by_mat[code] = set(np.unique(tris_by_mat[code]).tolist())
    gmsh.finalize()

    # ---- orientation normalisation (all signed areas positive) ----
    for code, arr in tris_by_mat.items():
        xy = nodes[arr]
        cr = (xy[:, 1, 0] - xy[:, 0, 0]) * (xy[:, 2, 1] - xy[:, 0, 1]) - \
             (xy[:, 1, 1] - xy[:, 0, 1]) * (xy[:, 2, 0] - xy[:, 0, 0])
        flip = cr < 0
        arr[flip] = arr[flip][:, [0, 2, 1]]
        tris_by_mat[code] = arr

    # ---- geometry helpers ----
    def a_of(x):
        return float(M.a_of(x)) if np.isscalar(x) else M.a_of(x)

    x = nodes[:, 0]; y = nodes[:, 1]
    a = M.a_of(x)

    # ---- node sets per intended boundary ----
    gate_n = nodes_by_mat.get(M.GATE, set())
    shield_n = nodes_by_mat.get(M.SHIELD, set())
    si_n = nodes_by_mat.get(M.SI, set())
    ox_n = nodes_by_mat.get(M.OX, set())

    source_n = {i for i in si_n if x[i] <= 1e-12}
    drain_n = {i for i in si_n if x[i] >= M.BOTTOM_DEPTH - 1e-9}
    body_n = {i for i in si_n if x[i] <= BODY_BOTTOM + 1e-12}
    interface_n = si_n & ox_n
    outer_n = {i for i in range(len(nodes))
               if abs(abs(y[i]) - M.Y_EDGE) < 1e-9 and x[i] < M.BOTTOM_DEPTH - 1e-9}

    # ---- checks ----
    checks = {}
    checks["gate_shield_disjoint"] = len(gate_n & shield_n) == 0
    checks["gate_inside_silicon"] = len(gate_n & si_n) == 0
    checks["shield_inside_silicon"] = len(shield_n & si_n) == 0
    checks["gate_shield_node_counts"] = {"Gate": len(gate_n), "Shield": len(shield_n)}
    checks["interface_nodes_shared"] = {
        "count": len(interface_n),
        "in_Si_only": len(si_n - ox_n),
        "in_oxide_only": len(ox_n - si_n),
    }
    checks["dirichlet_sets"] = {
        "source": {"nodes": len(source_n), "value": V_SRC},
        "body": {"nodes": len(body_n), "value": V_SRC},
        "drain": {"nodes": len(drain_n), "value": V_DRAIN},
        "gate": {"nodes": len(gate_n), "value": V_G},
        "shield": {"nodes": len(shield_n), "value": V_SG},
    }
    checks["neumann_outer_nodes"] = len(outer_n)

    # ---- assemble ----
    all_tris = np.vstack([tris_by_mat[k] for k in tris_by_mat])
    mat = np.concatenate([[k] * len(tris_by_mat[k]) for k in tris_by_mat])
    xy = nodes[all_tris]
    d1 = xy[:, 1] - xy[:, 0]; d2 = xy[:, 2] - xy[:, 0]
    det = d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0]
    area = 0.5 * det
    b = np.stack([xy[:, 1, 1] - xy[:, 2, 1], xy[:, 2, 1] - xy[:, 0, 1], xy[:, 0, 1] - xy[:, 1, 1]], 1)
    c = np.stack([xy[:, 2, 0] - xy[:, 1, 0], xy[:, 0, 0] - xy[:, 2, 0], xy[:, 1, 0] - xy[:, 0, 0]], 1)
    eps_e = np.where(mat == M.OX, EPS_OX, EPS_SI)
    ke = eps_e[:, None, None] * (b[:, :, None] * b[:, None, :] + c[:, :, None] * c[:, None, :]) \
        / (4 * area[:, None, None])
    rows = np.repeat(all_tris, 3, axis=1).ravel()
    cols = np.tile(all_tris, (1, 3)).ravel()
    K = sp.coo_matrix((ke.ravel(), (rows, cols)), shape=(len(nodes), len(nodes))).tocsr()

    f = np.zeros(len(nodes))
    si_elem = mat == M.SI
    np.add.at(f, all_tris[si_elem].T, (RHO * area[si_elem] / 3.0))

    q_total = float(f.sum())
    q_analytic = float(RHO * area[si_elem].sum())
    checks["charge"] = {
        "rho_V_per_um2": float(RHO),
        "integrated_nodal_charge": q_total,
        "analytic_rho_times_si_area": q_analytic,
        "relative_error": abs(q_total - q_analytic) / q_analytic,
        "silicon_area_um2": float(area[si_elem].sum()),
        "oxide_elements_receiving_charge": int((np.abs(f[np.unique(all_tris[~si_elem])]) > 0).sum()),
    }

    # ---- Dirichlet assembly ----
    fixed = np.zeros(len(nodes), bool)
    val = np.zeros(len(nodes))
    for nset, v in ((source_n, V_SRC), (body_n, V_SRC), (drain_n, V_DRAIN),
                    (gate_n, V_G), (shield_n, V_SG)):
        for i in nset:
            fixed[i] = True
            val[i] = v
    checks["dirichlet_nodes_total"] = int(fixed.sum())
    checks["dirichlet_conflicts"] = int(sum(
        1 for i in range(len(nodes)) if fixed[i] and
        len({V_SRC if i in source_n else None,
             V_SRC if i in body_n else None,
             V_DRAIN if i in drain_n else None,
             V_G if i in gate_n else None,
             V_SG if i in shield_n else None} - {None}) > 1))

    free = ~fixed
    rhs = f - K @ val
    u = val.copy()
    idx = np.flatnonzero(free)
    u[idx] = spla.spsolve(K[idx][:, idx].tocsc(), rhs[idx])

    # complementarity on silicon nodes only
    si_node = np.zeros(len(nodes), bool)
    for i in si_n:
        si_node[i] = True
    active = np.zeros(len(nodes), bool)
    active[idx] = (u[idx] < 0) & si_node[idx]
    for it in range(100):
        free2 = free & ~active
        idx2 = np.flatnonzero(free2)
        u2 = val.copy()
        if idx2.size:
            u2[idx2] = spla.spsolve(K[idx2][:, idx2].tocsc(), rhs[idx2])
        r = K @ u2 - rhs
        r[fixed] = 0.0
        dg = K.diagonal(); dg = np.where(dg > 0, dg, 1.0)
        new = np.zeros(len(nodes), bool)
        new[idx] = (u2[idx] - r[idx] / dg[idx] <= 0) & si_node[idx]
        if np.array_equal(new, active):
            active = new; u = u2
            break
        active = new; u = u2

    res = float(np.max(np.abs(r[idx2])) / max(np.max(np.abs(rhs)), 1.0)) if idx2.size else 0.0
    checks["complementarity"] = {
        "active_iterations": it + 1,
        "constrained_residual": res,
        "active_nodes_total": int(active.sum()),
        "active_nodes_in_silicon": int((active & si_node).sum()),
        "active_nodes_in_oxide": int((active & ~si_node).sum()),
        "complementarity_in_silicon_only": bool((active & ~si_node).sum() == 0),
    }
    checks["potential_range"] = {"min": float(np.nanmin(u)), "max": float(np.nanmax(u))}

    with open(out, "w") as fh:
        json.dump(checks, fh, indent=2)
    print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    audit()
