#!/usr/bin/env python3
"""Gmsh unstructured conforming mesh + P1 FEM Poisson with a silicon-only
complementarity constraint.  Vectorised assembly and KD-tree probes so the
three-resolution gate is tractable.

Physics: div(eps grad phi) = -rho in Si (rho = qND), 0 in oxide.
Silicon obstacle: u = V_DRAIN - phi >= 0 with -lap u + qND/eps >= 0.
Mesa capped at P/2 - t_ox so the oxide keeps full thickness through the tip.
"""
from __future__ import annotations
import argparse, json, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.spatial import cKDTree
import gmsh

PITCH = 1.35
A_TOP, A_BOT = 0.34, 0.425
TRENCH_DEPTH, TIP_START, OXIDE_T = 2.20, 2.05, 0.10
BOTTOM_DEPTH, BODY_BOTTOM = 2.42, 0.54
GATE_TOP, GATE_BOTTOM = 0.10, 0.54
SHIELD_TOP, SHIELD_BOTTOM = 0.69, 2.04
A_MAX = PITCH / 2 - OXIDE_T

ND = 6.9e16
EPS_SI, EPS_OX = 11.7, 3.9
EPS0_FCM = 8.8541878128e-14
Q = 1.602176634e-19
V_DRAIN = 34.274208536
RHO_SI = (Q * ND / EPS0_FCM) * 1e-8

SI, OX, METAL = 0, 1, 2


def half_width(x, r_tip):
    x = np.asarray(x, float)
    taper = A_TOP + (A_BOT - A_TOP) * np.clip(x / TIP_START, 0.0, 1.0)
    r = max(float(r_tip), 1e-9)
    u = np.clip((x - TIP_START) / r, 0.0, 1.0)
    arc = A_BOT + r - r * np.sqrt(np.maximum(1.0 - u * u, 0.0))
    xa, xb = TIP_START + r, TRENCH_DEPTH
    blend = np.clip((x - xa) / max(xb - xa, 1e-12), 0.0, 1.0)
    cap = A_BOT + r + (A_MAX - (A_BOT + r)) * blend ** 2
    return np.where(x <= TIP_START, taper, np.where(x <= xa, arc, np.minimum(cap, A_MAX)))


def build_mesh(h_tip, h_far, r_tip):
    t0 = time.time()
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 1)
    gmsh.model.add("trench")
    occ = gmsh.model.occ
    rect = occ.addRectangle(0.0, -PITCH / 2, 0.0, BOTTOM_DEPTH, PITCH)
    occ.synchronize()
    print(f"[mesh] rectangle {time.time()-t0:.1f}s", flush=True)

    xs = np.linspace(0.0, BOTTOM_DEPTH, 48)
    a = half_width(xs, r_tip)
    curves = []
    for sgn in (+1.0, -1.0):
        pts = [occ.addPoint(float(xs[i]), float(sgn * a[i]), 0.0, h_far)
               for i in range(len(xs))]
        curves.append(occ.addSpline(pts))
    occ.synchronize()
    for c in curves:
        gmsh.model.mesh.embed(1, [c], 2, rect)
    occ.synchronize()
    print(f"[mesh] curves embedded {time.time()-t0:.1f}s", flush=True)

    f = gmsh.model.mesh.field
    box = f.add("Box")
    f.setNumber(box, "VIn", h_tip); f.setNumber(box, "VOut", h_far)
    f.setNumber(box, "XMin", TIP_START - 0.15); f.setNumber(box, "XMax", TRENCH_DEPTH + 0.05)
    f.setNumber(box, "YMin", -A_MAX - 0.15); f.setNumber(box, "YMax", A_MAX + 0.15)
    f.setNumber(box, "Thickness", 0.25)
    f.setAsBackgroundMesh(box)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    gmsh.model.mesh.generate(2)
    print(f"[mesh] generated {time.time()-t0:.1f}s", flush=True)

    ntags, ncoords, _ = gmsh.model.mesh.getNodes()
    nodes = np.column_stack([ncoords[0::3], ncoords[1::3]])
    tag2i = {int(t): i for i, t in enumerate(ntags)}
    etypes, _, enodes = gmsh.model.mesh.getElements(dim=2)
    tris = np.vstack([np.asarray(en, dtype=int).reshape(-1, 3)
                      for et, en in zip(etypes, enodes) if et == 2])
    tris = np.vectorize(tag2i.get)(tris)
    gmsh.finalize()
    print(f"[mesh] {len(nodes)} nodes, {len(tris)} triangles, {time.time()-t0:.1f}s", flush=True)
    return nodes, tris


def materials(nodes, tris, r_tip):
    cx = nodes[tris, 0].mean(axis=1); cy = nodes[tris, 1].mean(axis=1)
    a = half_width(cx, r_tip)
    mat = np.where(np.abs(cy) < a, SI, OX)
    ay = np.abs(cy)
    metal = ((cx >= GATE_TOP) & (cx <= GATE_BOTTOM) & (ay >= a + OXIDE_T)) | \
            ((cx >= SHIELD_TOP) & (cx <= SHIELD_BOTTOM) & (ay >= a + OXIDE_T)) | \
            ((cx <= BODY_BOTTOM) & (np.abs(cy) < a))
    mat = np.where(metal, METAL, mat)
    return mat


def geom(nodes, tris):
    xy = nodes[tris]
    d = xy[:, 1] - xy[:, 0]; e = xy[:, 2] - xy[:, 0]
    det = d[:, 0] * e[:, 1] - d[:, 1] * e[:, 0]
    area = 0.5 * np.abs(det)
    b = np.stack([xy[:, 1, 1] - xy[:, 2, 1], xy[:, 2, 1] - xy[:, 0, 1], xy[:, 0, 1] - xy[:, 1, 1]], 1)
    c = np.stack([xy[:, 2, 0] - xy[:, 1, 0], xy[:, 0, 0] - xy[:, 2, 0], xy[:, 1, 0] - xy[:, 0, 0]], 1)
    return xy, det, area, b, c


def assemble(nodes, tris, mat):
    t0 = time.time()
    xy, det, area, b, c = geom(nodes, tris)
    eps = np.where(mat == OX, EPS_OX, EPS_SI)
    ke = eps[:, None, None] * (b[:, :, None] * b[:, None, :] + c[:, :, None] * c[:, None, :]) \
        / (4 * area[:, None, None])
    rows = np.repeat(tris, 3, axis=1).ravel()
    cols = np.tile(tris, (1, 3)).ravel()
    K = sp.coo_matrix((ke.ravel(), (rows, cols)), shape=(len(nodes), len(nodes))).tocsr()

    f = np.zeros(len(nodes))
    si = mat == SI
    contrib = (RHO_SI * area[si] / 3.0)
    np.add.at(f, tris[si].T, contrib)
    print(f"[assemble] K nnz={K.nnz} {time.time()-t0:.1f}s", flush=True)
    return K, f, area


def solve(K, f, fixed, val, si_any, max_iter=200):
    t0 = time.time()
    n = len(f)
    free = ~fixed
    rhs = f - K @ val
    u = val.copy()
    idx = np.flatnonzero(free)
    u[idx] = spla.spsolve(K[idx][:, idx].tocsc(), rhs[idx])
    active = np.zeros(n, bool)
    active[idx] = (u[idx] < 0) & si_any[idx]
    it = 0
    for it in range(max_iter):
        free2 = free & ~active
        idx2 = np.flatnonzero(free2)
        u2 = val.copy()
        if idx2.size:
            u2[idx2] = spla.spsolve(K[idx2][:, idx2].tocsc(), rhs[idx2])
        r = K @ u2 - rhs
        r[fixed] = 0.0
        dg = K.diagonal(); dg = np.where(dg > 0, dg, 1.0)
        new = np.zeros(n, bool)
        new[idx] = (u2[idx] - r[idx] / dg[idx] <= 0) & si_any[idx]
        if np.array_equal(new, active):
            active = new; u = u2
            break
        active = new; u = u2
    res = float(np.max(np.abs(r[idx2])) / max(np.max(np.abs(rhs)), 1.0)) if idx2.size else 0.0
    print(f"[solve] {it+1} active-set iterations, {time.time()-t0:.1f}s", flush=True)
    return u, active, it + 1, res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h-tip", type=float, default=0.02)
    ap.add_argument("--h-far", type=float, default=0.12)
    ap.add_argument("--r-tip", type=float, default=0.10)
    ap.add_argument("--prefix", default="fem")
    args = ap.parse_args()

    nodes, tris = build_mesh(args.h_tip, args.h_far, args.r_tip)
    mat = materials(nodes, tris, args.r_tip)
    K, f, area = assemble(nodes, tris, mat)

    x = nodes[:, 0]
    fixed_metal = np.zeros(len(nodes), bool)
    for e in np.flatnonzero(mat == METAL):
        fixed_metal[tris[e]] = True
    a = half_width(x, args.r_tip)
    src = (x <= 1e-12) & (np.abs(nodes[:, 1]) < a)
    drn = x >= BOTTOM_DEPTH - 1e-12
    fixed = fixed_metal | src | drn
    val = np.zeros(len(nodes)); val[drn] = V_DRAIN

    si_any = np.zeros(len(nodes), bool)
    for e in np.flatnonzero(mat == SI):
        si_any[tris[e]] = True

    u, active, iters, res = solve(K, f, fixed, val, si_any)
    phi = V_DRAIN - u

    # probes via KD-tree on element centroids
    cent = nodes[tris].mean(axis=1)
    tree = cKDTree(cent)
    a130 = float(half_width(1.30, args.r_tip))

    def field_at(px, py):
        for e in tree.query([px, py], k=24)[1]:
            xy = nodes[tris[e]]
            M = np.array([[xy[1, 0] - xy[0, 0], xy[2, 0] - xy[0, 0]],
                          [xy[1, 1] - xy[0, 1], xy[2, 1] - xy[0, 1]]])
            try:
                st = np.linalg.solve(M, np.array([px - xy[0, 0], py - xy[0, 1]]))
            except np.linalg.LinAlgError:
                continue
            if st[0] >= -1e-9 and st[1] >= -1e-9 and st[0] + st[1] <= 1 + 1e-9:
                det = np.linalg.det(M)
                b = np.array([xy[1, 1] - xy[2, 1], xy[2, 1] - xy[0, 1], xy[0, 1] - xy[1, 1]])
                c = np.array([xy[2, 0] - xy[1, 0], xy[0, 0] - xy[2, 0], xy[1, 0] - xy[0, 0]])
                gx = -(phi[tris[e]] @ b) / det; gy = -(phi[tris[e]] @ c) / det
                return float(np.hypot(gx, gy) * 0.01), float(phi[tris[e]] @ np.array(
                    [1 - st[0] - st[1], st[0], st[1]]))
        return float("nan"), float("nan")

    probes = {"Si_centre_1p30": (1.30, 0.0),
              "Si_inside_1p30": (1.30, a130 - 0.005),
              "Ox_inside_1p30": (1.30, a130 + 0.005)}
    out = {k: dict(zip(("E_MV_cm", "phi_V"), field_at(*p))) for k, p in probes.items()}

    # quality
    xy = nodes[tris]
    e_len = np.stack([np.linalg.norm(xy[:, 1] - xy[:, 0], axis=1),
                      np.linalg.norm(xy[:, 2] - xy[:, 1], axis=1),
                      np.linalg.norm(xy[:, 0] - xy[:, 2], axis=1)], 1)
    ang = []
    for k in range(3):
        v1 = xy[:, (k + 1) % 3] - xy[:, k]; v2 = xy[:, (k + 2) % 3] - xy[:, k]
        cosang = np.sum(v1 * v2, 1) / (np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1) + 1e-30)
        ang.append(np.degrees(np.arccos(np.clip(cosang, -1, 1))))
    ang = np.concatenate(ang)

    si_area = area[mat == SI].sum()
    neutral_area = area[(mat == SI) & active[tris[(mat == SI)]].any(axis=1)].sum()

    info = dict(
        geometry=dict(r_tip_um=args.r_tip, a_max_um=A_MAX, oxide_um=OXIDE_T),
        mesh=dict(vertices=int(len(nodes)), triangles=int(len(tris)),
                  min_angle_deg=float(ang.min()), mean_angle_deg=float(ang.mean()),
                  min_area_um2=float(area.min()),
                  h_tip=args.h_tip, h_far=args.h_far),
        solve=dict(iterations=iters, residual=res,
                   neutral_si_area_fraction=float(neutral_area / si_area)),
        probes=out,
        rho_si_V_per_um2=float(RHO_SI),
    )
    with open(f"{args.prefix}_check.json", "w") as fh:
        json.dump(info, fh, indent=2)
    print(json.dumps(info, indent=2))


if __name__ == "__main__":
    main()
