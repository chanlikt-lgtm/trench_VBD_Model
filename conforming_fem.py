#!/usr/bin/env python3
"""Gmsh unstructured conforming mesh + P1 finite-element Poisson with a
silicon-only complementarity (depletion) constraint.

Replaces the layered-band finite-volume version, which degenerated at the tip
(max non-orthogonality 79 deg, min area falling under refinement).

Physics kept identical in form:
    div(eps grad phi) = -rho,   rho = q ND in depleted silicon, 0 in oxide
    u = V_DRAIN - phi,  u >= 0,  -lap u + qND/eps >= 0,  u(...) = 0   in Si
    oxide: unconstrained, charge free, coupled through the same FE system.

Geometry fix required for a conforming mesh: the mesa half width is capped at
P/2 - t_ox so the oxide keeps its full thickness through the tip. The legacy
profile drove a(x) -> P/2, which makes the oxide zero-thickness there.
"""
from __future__ import annotations
import argparse
import json
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import gmsh

# ---------------------------------------------------------------- parameters
PITCH = 1.35
A_TOP = 0.34
A_BOT = 0.425
TRENCH_DEPTH = 2.20
TIP_START = 2.05
OXIDE_T = 0.10
BOTTOM_DEPTH = 2.42
BODY_BOTTOM = 0.54
GATE_TOP, GATE_BOTTOM = 0.10, 0.54
SHIELD_TOP, SHIELD_BOTTOM = 0.69, 2.04
A_MAX = PITCH / 2 - OXIDE_T          # 0.575: oxide never collapses

ND = 6.9e16
EPS_SI, EPS_OX = 11.7, 3.9
EPS0_FCM = 8.8541878128e-14
Q = 1.602176634e-19
V_DRAIN = 34.274208536
V_GROUND = 0.0
RHO_SI = (Q * ND / EPS0_FCM) * 1e-8    # V/um^2, verified 1.2486e3

SI, OX, METAL = 0, 1, 2


def half_width(x, r_tip):
    """Taper, circular tip arc of radius r_tip, smooth cap at A_MAX."""
    x = np.asarray(x, float)
    taper = A_TOP + (A_BOT - A_TOP) * np.clip(x / TIP_START, 0.0, 1.0)
    r = max(float(r_tip), 1e-9)
    u = np.clip((x - TIP_START) / r, 0.0, 1.0)
    arc = A_BOT + r - r * np.sqrt(np.maximum(1.0 - u * u, 0.0))
    xa, xb = TIP_START + r, TRENCH_DEPTH
    blend = np.clip((x - xa) / max(xb - xa, 1e-12), 0.0, 1.0)
    cap = A_BOT + r + (A_MAX - (A_BOT + r)) * blend ** 2
    return np.where(x <= TIP_START, taper, np.where(x <= xa, arc, np.minimum(cap, A_MAX)))


def build_mesh(h_tip, h_far, r_tip, verbose=False):
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 1 if verbose else 0)
    gmsh.model.add("trench")

    occ = gmsh.model.occ
    rect = occ.addRectangle(0.0, -PITCH / 2, 0.0, BOTTOM_DEPTH, PITCH)
    occ.synchronize()

    xs = np.linspace(0.0, BOTTOM_DEPTH, 400)
    a = half_width(xs, r_tip)

    curves = []
    for sgn in (+1.0, -1.0):
        pts = [occ.addPoint(float(xs[i]), float(sgn * a[i]), 0.0, h_far)
               for i in range(len(xs))]
        curves.append(occ.addSpline(pts))
    occ.synchronize()
    for c in curves:
        gmsh.model.mesh.embed(1, [c], 2, rect)

    # size field: fine at the tip, coarse far away
    f = gmsh.model.mesh.field
    box = f.add("Box")
    f.setNumber(box, "VIn", h_tip)
    f.setNumber(box, "VOut", h_far)
    f.setNumber(box, "XMin", TIP_START - 0.15)
    f.setNumber(box, "XMax", TRENCH_DEPTH + 0.05)
    f.setNumber(box, "YMin", -A_MAX - 0.15)
    f.setNumber(box, "YMax", A_MAX + 0.15)
    f.setNumber(box, "Thickness", 0.25)
    f.setAsBackgroundMesh(box)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.optimize("Netgen")

    ntags, ncoords, _ = gmsh.model.mesh.getNodes()
    nodes = np.column_stack([ncoords[0::3], ncoords[1::3]])
    tag2i = {int(t): i for i, t in enumerate(ntags)}

    etypes, etags, enodes = gmsh.model.mesh.getElements(dim=2)
    tris = []
    for et, en in zip(etypes, enodes):
        if et == 2:
            arr = np.asarray(en, dtype=int).reshape(-1, 3)
            tris.append(arr)
    tris = np.vstack(tris)
    tris = np.vectorize(tag2i.get)(tris)

    tri_nodes = np.asarray([gmsh.model.mesh.getNodes()[0]]).size  # noqa: F841
    gmsh.finalize()
    return nodes, tris


def element_material(nodes, tris, r_tip):
    """Material by element centroid: Si, oxide, or metal block."""
    cx = nodes[tris, 0].mean(axis=1)
    cy = nodes[tris, 1].mean(axis=1)
    a = half_width(cx, r_tip)
    mat = np.where(np.abs(cy) < a, SI, OX)
    ay = np.abs(cy)
    in_gate = (cx >= GATE_TOP) & (cx <= GATE_BOTTOM) & (ay >= a + OXIDE_T)
    in_shield = (cx >= SHIELD_TOP) & (cx <= SHIELD_BOTTOM) & (ay >= a + OXIDE_T)
    mat = np.where(in_gate | in_shield, METAL, mat)
    mat = np.where((cx <= BODY_BOTTOM) & (np.abs(cy) < a), METAL, mat)
    return mat


def assemble(nodes, tris, mat):
    n = len(nodes)
    rows, cols, vals = [], [], []
    for e in range(len(tris)):
        ids = tris[e]
        xy = nodes[ids]
        x1, y1 = xy[0]; x2, y2 = xy[1]; x3, y3 = xy[2]
        det = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
        area = 0.5 * abs(det)
        if area <= 0:
            continue
        b = np.array([y2 - y3, y3 - y1, y1 - y2])
        c = np.array([x3 - x2, x1 - x3, x2 - x1])
        eps = EPS_OX if mat[e] == OX else EPS_SI
        ke = eps * (np.outer(b, b) + np.outer(c, c)) / (4 * area)
        for i in range(3):
            for j in range(3):
                rows.append(ids[i]); cols.append(ids[j]); vals.append(ke[i, j])
    K = sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()

    # consistent load: rho * area / 3 per node of each Si element
    f = np.zeros(n)
    for e in range(len(tris)):
        if mat[e] != SI:
            continue
        ids = tris[e]
        xy = nodes[ids]
        det = (xy[1, 0] - xy[0, 0]) * (xy[2, 1] - xy[0, 1]) - \
              (xy[2, 0] - xy[0, 0]) * (xy[1, 1] - xy[0, 1])
        area = 0.5 * abs(det)
        f[ids] += RHO_SI * area / 3.0
    return K, f


def boundary_sets(nodes, tris, mat, r_tip):
    """Dirichlet node sets: source (x=0 Si), drain (x=BOTTOM), metal blocks."""
    n = len(nodes)
    is_metal = np.zeros(n, bool)
    for e in range(len(tris)):
        if mat[e] == METAL:
            is_metal[tris[e]] = True
    x = nodes[:, 0]
    a = half_width(x, r_tip)
    si = np.abs(nodes[:, 1]) < a
    src = (x <= 1e-12) & si
    drn = (x >= BOTTOM_DEPTH - 1e-12)
    fixed = is_metal | src | drn
    val = np.zeros(n)
    val[drn] = V_DRAIN
    return fixed, val


def solve_complementarity(K, f, fixed, val, si_nodes, max_iter=300):
    n = len(f)
    free = ~fixed
    u = np.zeros(n)
    rhs = f - K @ val
    Kf = K[free][:, free]
    rf = rhs[free]
    u_free = spla.spsolve(Kf.tocsc(), rf)
    u[:] = val
    u[free] = u_free

    active = np.zeros(n, bool)
    sf = si_nodes[free]
    active[free] = (u[free] < 0) & sf

    for it in range(max_iter):
        free2 = free & ~active
        u2 = np.zeros(n)
        u2[fixed] = val[fixed]
        idx = np.flatnonzero(free2)
        if idx.size:
            u2[idx] = spla.spsolve(K[idx][:, idx].tocsc(), rhs[idx])
        r = K @ u2 - rhs
        r[fixed] = 0.0
        diag = K.diagonal()
        diag = np.where(diag > 0, diag, 1.0)
        new = si_nodes & (u2 - r / diag <= 0) & free
        iters = it + 1
        if np.array_equal(new, active):
            break
        active = new
        u = u2
    else:
        iters = max_iter

    res = float(np.max(np.abs((K @ u2 - rhs)[free2])) / max(np.max(np.abs(rhs)), 1.0)) if free2.any() else 0.0
    return u2, active, iters, res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h-tip", type=float, default=0.005)
    ap.add_argument("--h-far", type=float, default=0.06)
    ap.add_argument("--r-tip", type=float, default=0.10)
    ap.add_argument("--prefix", default="fem")
    ap.add_argument("--save-tip-png", action="store_true")
    args = ap.parse_args()

    nodes, tris = build_mesh(args.h_tip, args.h_far, args.r_tip)
    mat = element_material(nodes, tris, args.r_tip)
    K, f = assemble(nodes, tris, mat)
    fixed, val = boundary_sets(nodes, tris, mat, args.r_tip)
    x = nodes[:, 0]
    a = half_width(x, args.r_tip)
    si_nodes = np.abs(nodes[:, 1]) < a
    matn = mat_of_nodes(nodes, tris, mat)
    si_nodes &= (matn == SI)
    # nodes on the Si/SiO2 interface belong to both materials; a node is silicon
    # if ANY incident element is silicon
    si_any = np.zeros(len(nodes), dtype=bool)
    for e in range(len(tris)):
        if mat[e] == SI:
            si_any[tris[e]] = True
    si_nodes = si_any

    u, active, iters, res = solve_complementarity(K, f, fixed, val, si_nodes)
    phi = V_DRAIN - u

    # mesh quality
    ang = []
    for e in range(len(tris)):
        xy = nodes[tris[e]]
        for k in range(3):
            v1 = xy[(k + 1) % 3] - xy[k]
            v2 = xy[(k + 2) % 3] - xy[k]
            c = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-30)
            ang.append(np.degrees(np.arccos(np.clip(c, -1, 1))))
    ang = np.array(ang)
    area = np.array([0.5 * abs(np.linalg.det(np.vstack([nodes[tris[e]][1] - nodes[tris[e]][0],
                                                        nodes[tris[e]][2] - nodes[tris[e]][0]])))
                     for e in range(len(tris))])

    # probes at fixed physical points
    def field_at(px, py):
        # local P1 gradient from the containing element
        for e in range(len(tris)):
            xy = nodes[tris[e]]
            v0 = xy[1] - xy[0]; v1 = xy[2] - xy[0]
            M = np.array([[v0[0], v1[0]], [v0[1], v1[1]]])
            try:
                st = np.linalg.solve(M, np.array([px - xy[0, 0], py - xy[0, 1]]))
            except np.linalg.LinAlgError:
                continue
            if st[0] >= -1e-9 and st[1] >= -1e-9 and st[0] + st[1] <= 1 + 1e-9:
                b = np.array([xy[1, 1] - xy[2, 1], xy[2, 1] - xy[0, 1], xy[0, 1] - xy[1, 1]])
                c = np.array([xy[2, 0] - xy[1, 0], xy[0, 0] - xy[2, 0], xy[1, 0] - xy[0, 0]])
                area2 = np.linalg.det(M)
                g = -(phi[tris[e]] @ b) / area2, -(phi[tris[e]] @ c) / area2
                return float(np.hypot(*g)) * 0.01, float(phi[tris[e]] @
                       np.array([1 - st[0] - st[1], st[0], st[1]]))
        return float("nan"), float("nan")

    a130 = float(half_width(1.30, args.r_tip))
    probes = {
        "Si_centre_1p30": (1.30, 0.0),
        "Si_inside_1p30": (1.30, a130 - 0.005),
        "Ox_inside_1p30": (1.30, a130 + 0.005),
    }
    probe_out = {}
    for k, (px, py) in probes.items():
        e, p = field_at(px, py)
        probe_out[k] = dict(E_MV_cm=e, phi_V=p)

    si_area = area[mat == SI].sum()
    neutral_area = area[(mat == SI)][active[tris[(mat == SI)]] .any(axis=1)].sum()
    front = np.nan
    j = np.argmin(np.abs(x - 1.30))
    ys = nodes[:, 1]
    sel = (np.abs(x - 1.30) < 0.01)
    if sel.any():
        front = float(np.max(np.abs(ys[sel][active[sel]]))) if active[sel].any() else 0.0

    info = dict(
        stage="Gmsh unstructured P1-FEM, silicon complementarity",
        geometry=dict(r_tip_um=args.r_tip, a_max_um=A_MAX, oxide_um=OXIDE_T,
                      note="a(x) capped at P/2 - t_ox so the oxide keeps 0.10 um at the tip"),
        mesh=dict(vertices=int(len(nodes)), triangles=int(len(tris)),
                  min_angle_deg=float(ang.min()), mean_angle_deg=float(ang.mean()),
                  min_area_um2=float(area.min()), max_area_um2=float(area.max())),
        solve=dict(obstacle_iterations=iters, residual=res,
                   neutral_si_area_fraction=float(neutral_area / si_area),
                   depletion_front_at_1p30_um=front),
        probes=probe_out,
        rho_si_V_per_um2=float(RHO_SI),
    )
    with open(f"{args.prefix}_check.json", "w") as fh:
        json.dump(info, fh, indent=2)
    print(json.dumps(info, indent=2))


def mat_of_nodes(nodes, tris, mat):
    out = np.full(len(nodes), -1, dtype=int)
    for e in range(len(tris)):
        out[tris[e]] = mat[e]
    return out


if __name__ == "__main__":
    main()
