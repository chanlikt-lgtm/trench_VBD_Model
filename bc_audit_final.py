#!/usr/bin/env python3
"""Stage 2B.1 (repeat) — compact the mesh, verify orphans, audit assembly, and
test the LINEAR two-material Poisson operator before enabling the obstacle.

Order of operations, deliberately:
  1. identify unused 2-D nodes and check none of them lies on a Dirichlet surface
  2. compact nodes/triangles, remap everything by index
  3. require: unused 2D nodes 0, zero free diagonals 0, non-finite entries 0,
     unanchored components 0, G/SG shared nodes 0, charge error ~0
  4. solve the LINEAR system only, report R_free
  5. only then run the silicon complementarity and report primal / dual /
     complementarity residuals separately
"""
from __future__ import annotations
import json
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.sparse.csgraph import connected_components

BODY_BOTTOM = 0.54
V_G, V_SG, V_DRAIN, V_SRC = 0.0, 0.0, 34.274208536, 0.0
ND = 6.9e16
EPS0_FCM = 8.8541878128e-14
Q = 1.602176634e-19
EPS_SI, EPS_OX = 11.7, 3.9
RHO = (Q * ND / EPS0_FCM) * 1e-8
A_TOP, A_BOT = 0.34, 0.425
TIP_START, R_TIP, X_ARC_END = 2.05, 0.10, 2.15
A_ARC_END = A_BOT + R_TIP
A_MAX = 0.575
BOTTOM_DEPTH, Y_EDGE = 2.42, 0.675
GATE_TOP, GATE_BOTTOM = 0.10, 0.54
SHIELD_TOP, SHIELD_BOTTOM = 0.69, 2.04

SI, OX, GATE, SHIELD = 0, 1, 2, 3


def a_of(x):
    x = np.asarray(x, float)
    taper = A_TOP + (A_BOT - A_TOP) * np.clip(x / TIP_START, 0, 1)
    u = np.clip((x - TIP_START) / R_TIP, 0, 1)
    arc = A_BOT + R_TIP - R_TIP * np.sqrt(np.maximum(1 - u * u, 0))
    t = np.clip((x - X_ARC_END) / (BOTTOM_DEPTH - X_ARC_END), 0, 1)
    return np.where(x <= TIP_START, taper,
                    np.where(x <= X_ARC_END, arc, A_ARC_END + (A_MAX - A_ARC_END) * t))


def assemble(nodes, tris, mats):
    xy = nodes[tris]
    d1 = xy[:, 1] - xy[:, 0]; d2 = xy[:, 2] - xy[:, 0]
    det = d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0]
    area = 0.5 * det
    b = np.stack([xy[:, 1, 1] - xy[:, 2, 1], xy[:, 2, 1] - xy[:, 0, 1],
                  xy[:, 0, 1] - xy[:, 1, 1]], 1)
    c = np.stack([xy[:, 2, 0] - xy[:, 1, 0], xy[:, 0, 0] - xy[:, 2, 0],
                  xy[:, 1, 0] - xy[:, 0, 0]], 1)
    eps = np.where(mats == OX, EPS_OX, EPS_SI)
    ke = eps[:, None, None] * (b[:, :, None] * b[:, None, :] + c[:, :, None] * c[:, None, :]) \
        / (4 * area[:, None, None])
    K = sp.coo_matrix((ke.ravel(),
                       (np.repeat(tris, 3, axis=1).ravel(), np.tile(tris, (1, 3)).ravel())),
                      shape=(len(nodes), len(nodes))).tocsr()
    f = np.zeros(len(nodes))
    si = mats == SI
    contrib = RHO * area[si] / 3.0
    np.add.at(f, tris[si].T, np.broadcast_to(contrib, (3, contrib.size)))
    return K, f, area


def main(npz="mesh2b_mesh.npz", out="bc_audit_final.json"):
    d = np.load(npz)
    nodes_raw = d["nodes"]
    parts, mlist = [], []
    for code, key in ((SI, "Si"), (OX, "SiO2"), (GATE, "Gate"), (SHIELD, "Shield")):
        if key in d:
            arr = np.array(d[key], dtype=int)
            xy = nodes_raw[arr]
            cr = (xy[:, 1, 0] - xy[:, 0, 0]) * (xy[:, 2, 1] - xy[:, 0, 1]) - \
                 (xy[:, 1, 1] - xy[:, 0, 1]) * (xy[:, 2, 0] - xy[:, 0, 0])
            flip = cr < 0
            arr[flip] = arr[flip][:, [0, 2, 1]]
            parts.append(arr); mlist.append(np.full(len(arr), code))
    tris_raw = np.vstack(parts); mats = np.concatenate(mlist)

    # ---------- 1. orphans, checked BEFORE removal ----------
    used, inverse = np.unique(tris_raw.ravel(), return_inverse=True)
    orphans = np.setdiff1d(np.arange(len(nodes_raw)), used)
    ox_orph = nodes_raw[orphans] if len(orphans) else np.zeros((0, 2))
    on_dirichlet = []
    for k, (px, py) in enumerate(ox_orph):
        on_src = px <= 1e-12
        on_drn = px >= BOTTOM_DEPTH - 1e-9
        in_gate = (GATE_TOP <= px <= GATE_BOTTOM) and abs(py) >= float(a_of(px)) + 0.10 - 1e-6
        in_shield = (SHIELD_TOP <= px <= SHIELD_BOTTOM) and abs(py) >= float(a_of(px)) + 0.10 - 1e-6
        if on_src or on_drn or in_gate or in_shield:
            on_dirichlet.append(dict(index=int(orphans[k]), x=float(px), y=float(py),
                                     source=bool(on_src), drain=bool(on_drn),
                                     gate=bool(in_gate), shield=bool(in_shield)))

    nodes = nodes_raw[used]
    tris = inverse.reshape(tris_raw.shape)

    old_to_new = np.full(len(nodes_raw), -1, dtype=int)
    old_to_new[used] = np.arange(len(used))

    # ---------- 2/3. compact mesh state + boundary sets ----------
    nodes_of = {k: set() for k in (SI, OX, GATE, SHIELD)}
    for code in (SI, OX, GATE, SHIELD):
        sel = mats == code
        if sel.any():
            nodes_of[code] = set(np.unique(tris[sel]).tolist())

    x, y = nodes[:, 0], nodes[:, 1]
    si_n, ox_n = nodes_of[SI], nodes_of[OX]
    gate_n, shield_n = nodes_of[GATE], nodes_of[SHIELD]
    source_n = {i for i in si_n if x[i] <= 1e-12}
    drain_n = {i for i in si_n if x[i] >= BOTTOM_DEPTH - 1e-9}
    body_n = {i for i in si_n if x[i] <= BODY_BOTTOM + 1e-12}
    interface_n = si_n & ox_n
    outer_n = {i for i in range(len(nodes))
               if abs(abs(y[i]) - Y_EDGE) < 1e-9 and x[i] < BOTTOM_DEPTH - 1e-9}

    K, f, area = assemble(nodes, tris, mats)

    fixed = np.zeros(len(nodes), bool); val = np.zeros(len(nodes))
    for nset, v in ((source_n, V_SRC), (body_n, V_SRC), (drain_n, V_DRAIN),
                    (gate_n, V_G), (shield_n, V_SG)):
        for i in nset:
            fixed[i] = True; val[i] = v
    free = ~fixed
    idx = np.flatnonzero(free)
    dg = K.diagonal()

    # connected components over the whole mesh, then anchor test
    rows, cols = [], []
    for e in range(len(tris)):
        a_, b_, c_ = tris[e]
        rows += [a_, b_, b_, c_, c_, a_]
        cols += [b_, a_, c_, b_, a_, c_]
    A = sp.coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(nodes), len(nodes)))
    ncomp, labels = connected_components(A.tocsr(), directed=False)
    unanchored = 0
    for c in range(ncomp):
        sel = labels == c
        if not (sel & fixed).any():
            unanchored += 1

    engine = {
        "orphan_nodes_before": int(len(orphans)),
        "orphans_on_dirichlet_surface": on_dirichlet,
        "nodes_after_compaction": int(len(nodes)),
        "triangles": int(len(tris)),
        "unused_2d_nodes_after": int(len(np.setdiff1d(np.arange(len(nodes)), np.unique(tris)))),
        "free_nodes": int(idx.size),
        "free_nodes_zero_diagonal": int((dg[idx] <= 0).sum()),
        "min_free_diagonal": float(dg[idx].min()) if idx.size else None,
        "nonfinite_matrix_entries": int((~np.isfinite(K.data)).sum() + (~np.isfinite(f)).sum()),
        "connected_components": int(ncomp),
        "unanchored_components": int(unanchored),
        "gate_shield_shared_nodes": int(len(gate_n & shield_n)),
        "charge_relative_error": float(abs(f.sum() - RHO * area[mats == SI].sum())
                                       / (RHO * area[mats == SI].sum())),
    }
    engine["pass"] = bool(
        engine["free_nodes_zero_diagonal"] == 0 and
        engine["nonfinite_matrix_entries"] == 0 and
        engine["unanchored_components"] == 0 and
        engine["gate_shield_shared_nodes"] == 0 and
        engine["unused_2d_nodes_after"] == 0)

    info = dict(stage="2B.1 repeat: compacted mesh + linear operator test",
                engine=engine)

    if not engine["pass"]:
        info["status"] = "BLOCKED at assembly gate; no solve attempted"
        with open(out, "w") as fh:
            json.dump(info, fh, indent=2)
        print(json.dumps(info, indent=2))
        return

    # ---------- 4. linear operator test only ----------
    fixed_idx = np.flatnonzero(fixed)
    KFF = K[idx][:, idx].tocsc()
    KFD = K[idx][:, fixed_idx]
    phi_D = val[fixed_idx]
    fF = f[idx]
    phi_F = spla.spsolve(KFF, fF - KFD @ phi_D)
    num = float(np.linalg.norm(KFF @ phi_F + KFD @ phi_D - fF))
    den = float(np.linalg.norm(KFF @ phi_F) + np.linalg.norm(KFD @ phi_D) + np.linalg.norm(fF))
    info["linear_operator"] = {
        "R_free": num / den if den else float("nan"),
        "finite": bool(np.all(np.isfinite(phi_F))),
        "phi_min": float(np.min(phi_F)), "phi_max": float(np.max(phi_F)),
    }

    # ---------- 5. silicon complementarity ----------
    u = val.copy(); u[idx] = phi_F
    si_node = np.zeros(len(nodes), bool)
    for i in si_n:
        si_node[i] = True
    active = np.zeros(len(nodes), bool)
    active[idx] = (u[idx] < 0) & si_node[idx]
    it = 0
    for it in range(200):
        free2 = ~fixed & ~active
        idx2 = np.flatnonzero(free2)
        u2 = val.copy()
        if idx2.size:
            u2[idx2] = spla.spsolve(K[idx2][:, idx2].tocsc(), f[idx2] - K[idx2][:, fixed_idx] @ phi_D)
        r = K @ u2 - f
        r[fixed] = 0.0
        dgs = np.where(dg > 0, dg, 1.0)
        new = np.zeros(len(nodes), bool)
        new[idx] = (u2[idx] - r[idx] / dgs[idx] <= 0) & si_node[idx]
        if np.array_equal(new, active):
            active = new; u = u2; break
        active = new; u = u2

    r = K @ u - f; r[fixed] = 0.0
    up = np.maximum(-u, 0.0)
    act = active & free
    info["complementarity"] = {
        "iterations": it + 1,
        "primal_max_violation": float(up.max()),
        "dual_max_violation": float(max(0.0, -(r[act]).max()) if act.any() else 0.0),
        "complementarity_max": float(np.max(np.abs(u[act] * r[act]))) if act.any() else 0.0,
        "free_residual_rel": float(np.linalg.norm(r[free]) / (np.linalg.norm(f) + 1e-30)),
        "active_nodes": int(active.sum()),
        "active_in_oxide": int((active & ~si_node).sum()),
        "finite": bool(np.all(np.isfinite(u))),
    }
    info["potential"] = {"min": float(np.nanmin(u)), "max": float(np.nanmax(u))}

    with open(out, "w") as fh:
        json.dump(info, fh, indent=2)
    print(json.dumps(info, indent=2))


if __name__ == "__main__":
    main()
