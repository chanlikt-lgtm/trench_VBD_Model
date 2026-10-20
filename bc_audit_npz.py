#!/usr/bin/env python3
"""Stage 2B.1 — assembled boundary-condition and charge audit.

Reads the validated Stage 2B mesh from `mesh2b_mesh.npz` (written by
mesh_stage2b.py) so that no Gmsh call happens in this process. Audits what is
actually assembled into the P1 system, not just what the CAD tags claim.
"""
from __future__ import annotations
import json
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

BODY_BOTTOM = 0.54
V_G, V_SG, V_DRAIN, V_SRC = 0.0, 0.0, 34.274208536, 0.0
ND = 6.9e16
EPS0_FCM = 8.8541878128e-14
Q = 1.602176634e-19
EPS_SI, EPS_OX = 11.7, 3.9
RHO = (Q * ND / EPS0_FCM) * 1e-8
A_TOP, A_BOT = 0.34, 0.425
TIP_START, R_TIP, X_ARC_END = 2.05, 0.10, 2.15
A_MAX = 0.575
BOTTOM_DEPTH, Y_EDGE = 2.42, 0.675

SI, OX, GATE, SHIELD = 0, 1, 2, 3
NAMES = {SI: "Si", OX: "SiO2", GATE: "Gate", SHIELD: "Shield"}


def a_of(x):
    x = np.asarray(x, float)
    taper = A_TOP + (A_BOT - A_TOP) * np.clip(x / TIP_START, 0, 1)
    u = np.clip((x - TIP_START) / R_TIP, 0, 1)
    arc = A_BOT + R_TIP - R_TIP * np.sqrt(np.maximum(1 - u * u, 0))
    t = np.clip((x - X_ARC_END) / (BOTTOM_DEPTH - X_ARC_END), 0, 1)
    return np.where(x <= TIP_START, taper, np.where(x <= X_ARC_END, arc,
                    A_ARC_END + (A_MAX - A_ARC_END) * t))


A_ARC_END = A_BOT + R_TIP


def main(npz="mesh2b_mesh.npz", out="bc_audit.json"):
    d = np.load(npz)
    nodes = d["nodes"]
    tris, mats = [], []
    for code, key in ((SI, "Si"), (OX, "SiO2"), (GATE, "Gate"), (SHIELD, "Shield")):
        if key not in d:
            continue
        arr = np.array(d[key], dtype=int)
        xy = nodes[arr]
        cr = (xy[:, 1, 0] - xy[:, 0, 0]) * (xy[:, 2, 1] - xy[:, 0, 1]) - \
             (xy[:, 1, 1] - xy[:, 0, 1]) * (xy[:, 2, 0] - xy[:, 0, 0])
        flip = cr < 0
        arr[flip] = arr[flip][:, [0, 2, 1]]
        tris.append(arr); mats.append(np.full(len(arr), code))
    tris = np.vstack(tris); mats = np.concatenate(mats)

    nodes_of = {k: set() for k in NAMES}
    for code, arr in zip(mats, []):
        pass
    for k in NAMES:
        sel = mats == k
        if sel.any():
            nodes_of[k] = set(np.unique(tris[sel]).tolist())

    x, y = nodes[:, 0], nodes[:, 1]
    a = a_of(x)

    gate_n, shield_n = nodes_of[GATE], nodes_of[SHIELD]
    si_n, ox_n = nodes_of[SI], nodes_of[OX]
    source_n = {i for i in si_n if x[i] <= 1e-12}
    drain_n = {i for i in si_n if x[i] >= BOTTOM_DEPTH - 1e-9}
    body_n = {i for i in si_n if x[i] <= BODY_BOTTOM + 1e-12}
    interface_n = si_n & ox_n
    outer_n = {i for i in range(len(nodes))
               if abs(abs(y[i]) - Y_EDGE) < 1e-9 and x[i] < BOTTOM_DEPTH - 1e-9}

    info = {"stage": "2B.1 assembled BC + charge audit"}
    info["groups"] = {
        "Si": {"faces": 1, "nodes": len(si_n), "condition": "semiconductor (unknown)"},
        "SiO2": {"faces": 2, "nodes": len(ox_n), "condition": "dielectric (unknown)"},
        "Gate": {"faces": 2, "nodes": len(gate_n), "condition": f"Dirichlet V_G={V_G}"},
        "Shield": {"faces": 2, "nodes": len(shield_n), "condition": f"Dirichlet V_SG={V_SG}"},
    }
    info["electrode_isolation"] = {
        "gate_shield_shared_nodes": len(gate_n & shield_n),
        "gate_oxide_shared_nodes": len(gate_n & ox_n),
        "shield_oxide_shared_nodes": len(shield_n & ox_n),
        "gate_silicon_shared_nodes": len(gate_n & si_n),
        "shield_silicon_shared_nodes": len(shield_n & si_n),
        "electrically_separate": len(gate_n & shield_n) == 0,
    }
    info["dirichlet"] = {
        "source": {"nodes": len(source_n), "V": V_SRC},
        "body": {"nodes": len(body_n), "V": V_SRC},
        "drain": {"nodes": len(drain_n), "V": V_DRAIN},
        "gate": {"nodes": len(gate_n), "V": V_G},
        "shield": {"nodes": len(shield_n), "V": V_SG},
    }
    info["interface"] = {
        "shared_nodes": len(interface_n),
        "silicon_only_nodes": len(si_n - ox_n),
        "oxide_only_nodes": len(ox_n - si_n),
        "conforming": len(interface_n) > 0,
    }
    info["neumann"] = {"outer_pitch_edge_nodes": len(outer_n)}

    # ---------------- assemble ----------------
    xy = nodes[tris]
    d1 = xy[:, 1] - xy[:, 0]; d2 = xy[:, 2] - xy[:, 0]
    det = d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0]
    area = 0.5 * det
    b = np.stack([xy[:, 1, 1] - xy[:, 2, 1], xy[:, 2, 1] - xy[:, 0, 1], xy[:, 0, 1] - xy[:, 1, 1]], 1)
    c = np.stack([xy[:, 2, 0] - xy[:, 1, 0], xy[:, 0, 0] - xy[:, 2, 0], xy[:, 1, 0] - xy[:, 0, 0]], 1)
    eps_e = np.where(mats == OX, EPS_OX, EPS_SI)
    ke = eps_e[:, None, None] * (b[:, :, None] * b[:, None, :] + c[:, :, None] * c[:, None, :]) \
        / (4 * area[:, None, None])
    rows = np.repeat(tris, 3, axis=1).ravel()
    cols = np.tile(tris, (1, 3)).ravel()
    K = sp.coo_matrix((ke.ravel(), (rows, cols)), shape=(len(nodes), len(nodes))).tocsr()

    f = np.zeros(len(nodes))
    si_elem = mats == SI
    contrib = RHO * area[si_elem] / 3.0
    np.add.at(f, tris[si_elem].T, np.broadcast_to(contrib, (3, contrib.size)))
    info["charge"] = {
        "rho_V_per_um2": float(RHO),
        "silicon_area_um2": float(area[si_elem].sum()),
        "integrated_charge": float(f.sum()),
        "analytic_rho_area": float(RHO * area[si_elem].sum()),
        "relative_error": float(abs(f.sum() - RHO * area[si_elem].sum()) / (RHO * area[si_elem].sum())),
        "oxide_nodes_with_charge": int((np.abs(f[list(ox_n)]) > 0).sum()),
        "metal_nodes_with_charge": int((np.abs(f[list(gate_n | shield_n)]) > 0).sum()),
    }

    fixed = np.zeros(len(nodes), bool); val = np.zeros(len(nodes))
    for nset, v in ((source_n, V_SRC), (body_n, V_SRC), (drain_n, V_DRAIN),
                    (gate_n, V_G), (shield_n, V_SG)):
        for i in nset:
            fixed[i] = True; val[i] = v
    info["dirichlet"]["total_fixed_nodes"] = int(fixed.sum())
    info["dirichlet"]["free_nodes"] = int((~fixed).sum())

    rhs = f - K @ val
    u = val.copy()
    idx = np.flatnonzero(~fixed)
    dg_all = K.diagonal()
    info["linear_system"] = {
        "free_nodes": int(idx.size),
        "free_nodes_with_zero_diagonal": int((dg_all[idx] <= 0).sum()),
        "min_free_diagonal": float(dg_all[idx].min()) if idx.size else None,
        "matrix_nnz": int(K.nnz),
    }
    u[idx] = spla.spsolve(K[idx][:, idx].tocsc(), rhs[idx])

    si_node = np.zeros(len(nodes), bool)
    for i in si_n:
        si_node[i] = True
    active = np.zeros(len(nodes), bool)
    active[idx] = (u[idx] < 0) & si_node[idx]
    it = 0
    for it in range(150):
        free2 = ~fixed & ~active
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
            active = new; u = u2; break
        active = new; u = u2

    res = float(np.max(np.abs(r[idx2])) / max(np.max(np.abs(rhs)), 1.0)) if idx2.size else 0.0
    info["complementarity"] = {
        "iterations": it + 1,
        "scaled_residual": res,
        "active_nodes_total": int(active.sum()),
        "active_in_silicon": int((active & si_node).sum()),
        "active_in_oxide": int((active & ~si_node).sum()),
        "silicon_only": bool((active & ~si_node).sum() == 0),
    }
    info["potential"] = {"min": float(u.min()), "max": float(u.max()),
                         "at_source": float(u[list(source_n)].mean()) if source_n else None,
                         "at_drain": float(u[list(drain_n)].mean()) if drain_n else None,
                         "at_gate": float(u[list(gate_n)].mean()) if gate_n else None,
                         "at_shield": float(u[list(shield_n)].mean()) if shield_n else None}

    # ---------------- interface flux continuity (one-sided) ----------------
    edge = {}
    for e in range(len(tris)):
        for k in range(3):
            i, j = tris[e][k], tris[e][(k + 1) % 3]
            key = (min(i, j), max(i, j))
            edge.setdefault(key, []).append(e)
    grads = []
    for e in range(len(tris)):
        grads.append(((b[e] @ u[tris[e]]) / det[e], (c[e] @ u[tris[e]]) / det[e]))
    grads = np.array(grads)
    jumps = []
    for key, els in edge.items():
        if len(els) != 2:
            continue
        e0, e1 = els
        if {mats[e0], mats[e1]} != {SI, OX}:
            continue
        p0, p1 = nodes[key[0]], nodes[key[1]]
        t = p1 - p0; t /= np.linalg.norm(t)
        n = np.array([t[1], -t[0]])
        g0 = grads[e0]; g1 = grads[e1]
        f0 = EPS_SI * (g0[0] * n[0] + g0[1] * n[1])
        f1 = EPS_OX * (g1[0] * n[0] + g1[1] * n[1])
        jumps.append(abs(f0 - f1) / max(abs(f0), abs(f1), 1e-30))
    jumps = np.array(jumps) if jumps else np.array([np.nan])
    info["interface_flux"] = {
        "edges_checked": int(len(jumps)),
        "median_relative_jump": float(np.nanmedian(jumps)),
        "p95_relative_jump": float(np.nanpercentile(jumps, 95)),
        "max_relative_jump": float(np.nanmax(jumps)),
    }

    with open(out, "w") as fh:
        json.dump(info, fh, indent=2)
    print(json.dumps(info, indent=2))


if __name__ == "__main__":
    main()
