#!/usr/bin/env python3
"""Stage 2B.1 (corrected) — silicon obstacle posed on u = V_D - phi.

Transformation, before Dirichlet elimination:
    K phi = b_phi        ->        K u = f_u,   f_u = V_D * K @ 1 - b_phi
    u_D = V_D  on source / body / G / SG (all at 0 V)
    u_D = 0    on drain (at V_D)

Obstacle, on unconstrained silicon DOFs only (interface trace nodes included):
    u_i >= 0,  r_i = (K u - f_u)_i >= 0,  u_i r_i = 0
Oxide: unconstrained, r_i = 0.

Reports primal / dual / complementarity separately, the active silicon AREA,
and the net reaction on interface nodes (an unintended sheet charge shows up
there as a non-zero sum).
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
A_ARC_END = A_BOT + R_TIP
A_MAX, BOTTOM_DEPTH, Y_EDGE = 0.575, 2.42, 0.675

SI, OX, GATE, SHIELD = 0, 1, 2, 3
TOL = 1e-9


def load(npz):
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
    used, inv = np.unique(tris_raw.ravel(), return_inverse=True)
    return nodes_raw[used], inv.reshape(tris_raw.shape), mats


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
    fphi = np.zeros(len(nodes))
    si = mats == SI
    contrib = RHO * area[si] / 3.0
    np.add.at(fphi, tris[si].T, np.broadcast_to(contrib, (3, contrib.size)))
    return K, fphi, area


def main(npz="mesh2b_mesh.npz", out="bc_audit_u.json"):
    nodes, tris, mats = load(npz)
    K, fphi, area = assemble(nodes, tris, mats)
    n = len(nodes)

    nodes_of = {}
    for code in (SI, OX, GATE, SHIELD):
        sel = mats == code
        nodes_of[code] = set(np.unique(tris[sel]).tolist()) if sel.any() else set()

    x, y = nodes[:, 0], nodes[:, 1]
    si_n, ox_n = nodes_of[SI], nodes_of[OX]
    gate_n, shield_n = nodes_of[GATE], nodes_of[SHIELD]
    source_n = {i for i in si_n if x[i] <= 1e-12}
    drain_n = {i for i in si_n if x[i] >= BOTTOM_DEPTH - 1e-9}
    body_n = {i for i in si_n if x[i] <= BODY_BOTTOM + 1e-12}
    interface_n = si_n & ox_n

    # ---- transformation to u ----
    ones = np.ones(n)
    f_u = V_DRAIN * (K @ ones) - fphi

    fixed = np.zeros(n, bool)
    uD = np.zeros(n)
    for nset, phi_val in ((source_n, V_SRC), (body_n, V_SRC),
                          (gate_n, V_G), (shield_n, V_SG), (drain_n, V_DRAIN)):
        for i in nset:
            fixed[i] = True
            uD[i] = V_DRAIN - phi_val

    free = ~fixed
    idx = np.flatnonzero(free)
    KFF = K[idx][:, idx].tocsc()
    KFD = K[idx][:, fixed]
    uF_lin = spla.spsolve(KFF, f_u[idx] - KFD @ uD[fixed])
    lin = {
        "R_free_u": float(np.linalg.norm(KFF @ uF_lin + KFD @ uD[fixed] - f_u[idx]) /
                          (np.linalg.norm(KFF @ uF_lin) + np.linalg.norm(KFD @ uD[fixed]) +
                           np.linalg.norm(f_u[idx]))),
        "finite": bool(np.all(np.isfinite(uF_lin))),
    }

    u = uD.copy(); u[idx] = uF_lin
    si_node = np.zeros(n, bool)
    for i in si_n:
        si_node[i] = True

    # diagnostics the user asked for, on the LINEAR solution
    u_lin_si = u[idx][si_node[idx]]
    pre = {
        "silicon_overshoot_max_V": float(np.maximum(-u_lin_si, 0).max()) if u_lin_si.size else 0.0,
        "silicon_min_u_V": float(u_lin_si.min()) if u_lin_si.size else None,
        "overshoot_is_at_unconstrained_silicon": bool(
            u_lin_si.size and np.maximum(-u_lin_si, 0).max() > TOL),
    }

    # ---- active set on u, silicon only ----
    active = np.zeros(n, bool)
    active[idx] = (u[idx] < 0) & si_node[idx]
    it = 0
    for it in range(300):
        free2 = ~fixed & ~active
        idx2 = np.flatnonzero(free2)
        u2 = uD.copy()
        if idx2.size:
            u2[idx2] = spla.spsolve(K[idx2][:, idx2].tocsc(),
                                    f_u[idx2] - K[idx2][:, fixed] @ uD[fixed])
        r = K @ u2 - f_u
        r[fixed] = 0.0
        dg = K.diagonal(); dg = np.where(dg > 0, dg, 1.0)
        new = np.zeros(n, bool)
        new[idx] = (u2[idx] - r[idx] / dg[idx] <= 0) & si_node[idx]
        if np.array_equal(new, active):
            active = new; u = u2
            break
        active = new; u = u2

    r = K @ u - f_u
    r[fixed] = 0.0
    phi = V_DRAIN - u

    si_free = np.zeros(n, bool); si_free[list(si_n)] = True; si_free &= free
    ox_free = np.zeros(n, bool); ox_free[list(ox_n)] = True; ox_free &= free
    iface_free = np.zeros(n, bool)
    for i in interface_n:
        if free[i]:
            iface_free[i] = True

    # nodal silicon area, so the active region can be reported as an AREA
    nodal_area = np.zeros(n)
    si = mats == SI
    np.add.at(nodal_area, tris[si].T,
              np.broadcast_to(area[si] / 3.0, (3, int(si.sum()))))

    r_scale = float(np.max(np.abs(f_u))) or 1.0
    u_scale = float(np.max(np.abs(u))) or 1.0

    info = {
        "stage": "2B.1 corrected: obstacle on u = V_D - phi",
        "transform": {"f_u_scale": float(np.max(np.abs(f_u))),
                      "u_D_source_body_gate_shield": V_DRAIN - 0.0,
                      "u_D_drain": V_DRAIN - V_DRAIN},
        "linear_operator": lin,
        "linear_solution_diagnostics": pre,
        "complementarity": {
            "iterations": it + 1,
            "primal_min_u_silicon_V": float(u[si_free].min()) if si_free.any() else None,
            "primal_ok": bool(si_free.any() and u[si_free].min() >= -TOL),
            "dual_min_r_silicon": float(r[si_free].min()) if si_free.any() else None,
            "dual_ok": bool(si_free.any() and r[si_free].min() >= -TOL * r_scale),
            "complementarity_max_abs": float(np.max(np.abs(u[si_free] * r[si_free])))
            if si_free.any() else 0.0,
            "complementarity_normalised": float(np.max(np.abs(u[si_free] * r[si_free])) /
                                                (u_scale * r_scale)) if si_free.any() else 0.0,
            "oxide_equilibrium_rel": float(np.linalg.norm(r[ox_free]) / r_scale)
            if ox_free.any() else 0.0,
            "active_nodes": int(active.sum()),
            "active_in_oxide": int((active & ~si_node).sum()),
            "active_silicon_area_um2": float(nodal_area[active].sum()),
            "silicon_area_um2": float(nodal_area[si_node].sum()),
            "interface_reaction_sum": float(r[iface_free].sum()) if iface_free.any() else 0.0,
            "interface_reaction_normalised": float(abs(r[iface_free].sum()) / r_scale)
            if iface_free.any() else 0.0,
            "interface_free_nodes": int(iface_free.sum()),
        },
        "potential": {"min": float(phi.min()), "max": float(phi.max()),
                      "reconstruction_ok": bool(np.allclose(phi, V_DRAIN - u))},
    }
    with open(out, "w") as fh:
        json.dump(info, fh, indent=2)
    print(json.dumps(info, indent=2))


if __name__ == "__main__":
    main()
