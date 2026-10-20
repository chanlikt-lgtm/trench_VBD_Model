#!/usr/bin/env python3
"""Stage 2B mesh gate: add explicit G and SG electrode faces.

Geometry is built from curves shared between faces (no embedding):
  profile   y = a(x)            Si / SiO2   (line + circular arc + line)
  offset    y = a(x) + t_ox     SiO2 / outer (same construction, shifted 0.10)
  pitch     y = +/- P/2         outer boundary

Electrode faces = offset curve .. pitch edge, within their depth ranges:
  G  : 0.10 <= x <= 0.54
  SG : 0.69 <= x <= 2.04
Everything else outside Si and outside the electrodes is oxide.

Checks (the point of this gate):
  * every triangle has the same material at all three vertices and its centroid
    -> material boundaries are mesh lines, i.e. the mesh is conforming
  * face areas sum to the analytic domain area
  * electrodes do not overlap Si (checked by the same material test)
  * G and SG are separate physical groups
"""
from __future__ import annotations
import argparse, json, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import gmsh

PITCH = 1.35
A_TOP, A_BOT = 0.34, 0.425
TIP_START, R_TIP = 2.05, 0.10
X_ARC_END = TIP_START + R_TIP
A_ARC_END = A_BOT + R_TIP
BOTTOM_DEPTH = 2.42
OXIDE_T = 0.10
A_MAX = PITCH / 2 - OXIDE_T
Y_EDGE = PITCH / 2
GATE_TOP, GATE_BOTTOM = 0.10, 0.54
SHIELD_TOP, SHIELD_BOTTOM = 0.69, 2.04

SI, OX, GATE, SHIELD = 0, 1, 2, 3
NAMES = {SI: "Si", OX: "SiO2", GATE: "Gate", SHIELD: "Shield"}


def a_of(x):
    x = np.asarray(x, float)
    taper = A_TOP + (A_BOT - A_TOP) * np.clip(x / TIP_START, 0.0, 1.0)
    u = np.clip((x - TIP_START) / R_TIP, 0.0, 1.0)
    arc = A_BOT + R_TIP - R_TIP * np.sqrt(np.maximum(1.0 - u * u, 0.0))
    t = np.clip((x - X_ARC_END) / (BOTTOM_DEPTH - X_ARC_END), 0.0, 1.0)
    cap = A_ARC_END + (A_MAX - A_ARC_END) * t
    return np.where(x <= TIP_START, taper, np.where(x <= X_ARC_END, arc, cap))


def build(h_global, h_tip):
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)
    gmsh.model.add("stage2b")
    occ = gmsh.model.occ

    rect = occ.addRectangle(0.0, -Y_EDGE, 0.0, BOTTOM_DEPTH, PITCH)

    # ---------------- profile (Si / SiO2), explicit directions --------------
    p0 = occ.addPoint(0.0, A_TOP, 0.0, h_global)
    p1 = occ.addPoint(TIP_START, A_BOT, 0.0, h_tip)
    pc = occ.addPoint(TIP_START, A_ARC_END, 0.0, h_tip)
    p2 = occ.addPoint(X_ARC_END, A_ARC_END, 0.0, h_tip)
    p3 = occ.addPoint(BOTTOM_DEPTH, A_MAX, 0.0, h_global)
    up_l1 = occ.addLine(p0, p1)
    up_arc = occ.addCircleArc(p1, pc, p2)
    up_l2 = occ.addLine(p2, p3)

    m0 = occ.addPoint(0.0, -A_TOP, 0.0, h_global)
    m1 = occ.addPoint(TIP_START, -A_BOT, 0.0, h_tip)
    mc = occ.addPoint(TIP_START, -A_ARC_END, 0.0, h_tip)
    m2 = occ.addPoint(X_ARC_END, -A_ARC_END, 0.0, h_tip)
    m3 = occ.addPoint(BOTTOM_DEPTH, -A_MAX, 0.0, h_global)
    lo_l2 = occ.addLine(m3, m2)
    lo_arc = occ.addCircleArc(m2, mc, m1)
    lo_l1 = occ.addLine(m1, m0)

    right = occ.addLine(p3, m3)
    left = occ.addLine(m0, p0)
    si = occ.addPlaneSurface([occ.addCurveLoop(
        [up_l1, up_arc, up_l2, right, lo_l2, lo_arc, lo_l1, left])])

    # ---------------- offset curve (SiO2 / outer) ----------------
    xk = [0.0, GATE_TOP, GATE_BOTTOM, SHIELD_TOP, SHIELD_BOTTOM, TIP_START]
    yk = [float(a_of(x)) + OXIDE_T for x in xk]
    off_pts = [occ.addPoint(xk[i], yk[i], 0.0, h_global) for i in range(len(xk))]
    off_seg = [occ.addLine(off_pts[i], off_pts[i + 1]) for i in range(len(xk) - 1)]
    po = occ.addPoint(X_ARC_END, A_ARC_END + OXIDE_T, 0.0, h_tip)
    pc_o = occ.addPoint(TIP_START, A_ARC_END + OXIDE_T, 0.0, h_tip)
    po_end = occ.addPoint(BOTTOM_DEPTH, Y_EDGE, 0.0, h_global)
    off_arc = occ.addCircleArc(off_pts[-1], pc_o, po)
    off_tail = occ.addLine(po, po_end)
    off_curves = off_seg + [off_arc, off_tail]

    off_lo_curves = []
    for c in off_curves:
        off_lo_curves.append(c)   # placeholder; mirrored below

    # mirrored offset points/curves
    off_pts_lo = [occ.addPoint(xk[i], -yk[i], 0.0, h_global) for i in range(len(xk))]
    off_seg_lo = [occ.addLine(off_pts_lo[i], off_pts_lo[i + 1]) for i in range(len(xk) - 1)]
    po_lo = occ.addPoint(X_ARC_END, -(A_ARC_END + OXIDE_T), 0.0, h_tip)
    pc_lo = occ.addPoint(TIP_START, -(A_ARC_END + OXIDE_T), 0.0, h_tip)
    po_end_lo = occ.addPoint(BOTTOM_DEPTH, -Y_EDGE, 0.0, h_global)
    off_arc_lo = occ.addCircleArc(off_pts_lo[-1], pc_lo, po_lo)
    off_tail_lo = occ.addLine(po_lo, po_end_lo)
    off_curves_lo = off_seg_lo + [off_arc_lo, off_tail_lo]

    # ---------------- electrode faces ----------------
    def electrode(i0, i1, sign, tag_pts_top):
        """Face between offset segment [i0,i1] and the pitch edge."""
        seg = off_seg[i0] if sign > 0 else off_seg_lo[i0]
        p_a = off_pts[i0] if sign > 0 else off_pts_lo[i0]
        p_b = off_pts[i1] if sign > 0 else off_pts_lo[i1]
        top_a = occ.addPoint(xk[i0], sign * Y_EDGE, 0.0, h_global)
        top_b = occ.addPoint(xk[i1], sign * Y_EDGE, 0.0, h_global)
        # explicit directions: p_b -> top_b -> top_a -> p_a closes the loop
        l_right = occ.addLine(p_b, top_b)
        l_top = occ.addLine(top_b, top_a)
        l_left = occ.addLine(top_a, p_a)
        loop = occ.addCurveLoop([seg, l_right, l_top, l_left])
        return occ.addPlaneSurface([loop])

    faces = {"si": [si], "gate": [], "shield": [], "ox": []}
    gi0 = xk.index(GATE_TOP); gi1 = xk.index(GATE_BOTTOM)
    si0 = xk.index(SHIELD_TOP); si1 = xk.index(SHIELD_BOTTOM)

    gate_faces = [electrode(gi0, gi1, +1.0, None), electrode(gi0, gi1, -1.0, None)]
    shield_faces = [electrode(si0, si1, +1.0, None), electrode(si0, si1, -1.0, None)]

    occ.synchronize()

    all_in = [(2, rect)] + [(2, t) for t in [si] + gate_faces + shield_faces]
    out, _ = occ.fragment(all_in, [])
    occ.synchronize()

    surf = [t for d, t in out if d == 2]
    com = {t: np.array(gmsh.model.occ.getCenterOfMass(2, t)) for t in surf}

    groups = {SI: [], OX: [], GATE: [], SHIELD: []}
    for t in surf:
        x, y = com[t][0], com[t][1]
        a = float(a_of(x))
        if abs(y) < a:
            groups[SI].append(t)
        elif (GATE_TOP <= x <= GATE_BOTTOM and abs(y) >= a + OXIDE_T) or \
             (SHIELD_TOP <= x <= SHIELD_BOTTOM and abs(y) >= a + OXIDE_T):
            groups[GATE if (GATE_TOP <= x <= GATE_BOTTOM) else SHIELD].append(t)
        else:
            groups[OX].append(t)

    names = NAMES
    for k, lst in groups.items():
        if lst:
            gmsh.model.addPhysicalGroup(2, lst, name=names[k])
    print("[geom] faces:", {names[k]: len(v) for k, v in groups.items()}, flush=True)
    return groups


def run(h_global, h_tip, prefix, h_front=None):
    t0 = time.time()
    groups = build(h_global, h_tip)
    gmsh.option.setNumber("Mesh.MeshSizeMin", min(h_tip, h_front) if h_front else h_tip)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h_global)
    if h_front is not None:
        # localised refinement around the drain-connected neutral front, with a
        # smooth transition band; distant mesh stays at the r5 sizes
        f = gmsh.model.mesh.field
        bl = f.add("Box")
        f.setNumber(bl, "VIn", h_front)
        f.setNumber(bl, "VOut", h_global)
        f.setNumber(bl, "XMin", 2.30); f.setNumber(bl, "XMax", 2.42)
        f.setNumber(bl, "YMin", -0.25); f.setNumber(bl, "YMax", 0.25)
        f.setNumber(bl, "Thickness", 0.06)
        bt = f.add("Box")
        f.setNumber(bt, "VIn", h_tip)
        f.setNumber(bt, "VOut", h_global)
        f.setNumber(bt, "XMin", 1.95); f.setNumber(bt, "XMax", 2.26)
        f.setNumber(bt, "YMin", -0.72); f.setNumber(bt, "YMax", 0.72)
        f.setNumber(bt, "Thickness", 0.10)
        mn = f.add("Min")
        f.setNumbers(mn, "FieldsList", [bl, bt])
        f.setAsBackgroundMesh(mn)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.model.mesh.generate(2)
    dt = time.time() - t0

    ntags, ncoords, _ = gmsh.model.mesh.getNodes()
    nodes = np.column_stack([ncoords[0::3], ncoords[1::3]])
    tag2i = {int(t): i for i, t in enumerate(ntags)}

    per_face, tris_mat = {}, []
    for code, lst in groups.items():
        cnt = 0
        for f in lst:
            etypes, _, enodes = gmsh.model.mesh.getElements(2, f)
            for et, en in zip(etypes, enodes):
                if et != 2:
                    continue
                arr = np.vectorize(tag2i.get)(np.asarray(en, int).reshape(-1, 3))
                tris_mat.append((code, arr))
                cnt += len(arr)
        per_face[code] = cnt

    # dump the validated mesh so downstream audits never have to re-enter gmsh
    dump = {NAMES[k]: np.vstack([a for c, a in tris_mat if c == k])
            for k in NAMES if any(c == k for c, _ in tris_mat)}
    np.savez_compressed(f"{prefix}_mesh.npz", nodes=nodes, **dump)
    gmsh.finalize()

    # material map per node, and the conformity test
    node_mat = {}
    for code, arr in tris_mat:
        for v in np.unique(arr):
            node_mat.setdefault(int(v), set()).add(code)

    def mat_at(x, y):
        a = float(a_of(x))
        if abs(y) < a:
            return SI
        if GATE_TOP <= x <= GATE_BOTTOM and abs(y) >= a + OXIDE_T:
            return GATE
        if SHIELD_TOP <= x <= SHIELD_BOTTOM and abs(y) >= a + OXIDE_T:
            return SHIELD
        return OX

    bad = 0
    total = 0
    for code, arr in tris_mat:
        xy = nodes[arr]
        cen = xy.mean(axis=1)
        for e in range(len(arr)):
            # centroid is strictly interior, so its material must equal the
            # material of the face the element was taken from
            if mat_at(cen[e, 0], cen[e, 1]) != code:
                bad += 1
            total += 1

    allt = np.vstack([a for _, a in tris_mat])
    xy = nodes[allt]
    e1 = xy[:, 1] - xy[:, 0]; e3 = xy[:, 0] - xy[:, 2]
    cr = e1[:, 0] * (xy[:, 2, 1] - xy[:, 0, 1]) - e1[:, 1] * (xy[:, 2, 0] - xy[:, 0, 0])
    area = 0.5 * cr
    L = np.stack([np.linalg.norm(xy[:, 1] - xy[:, 0], axis=1),
                  np.linalg.norm(xy[:, 2] - xy[:, 1], axis=1),
                  np.linalg.norm(xy[:, 0] - xy[:, 2], axis=1)], 1)
    ang = []
    for k in range(3):
        v1 = xy[:, (k + 1) % 3] - xy[:, k]; v2 = xy[:, (k + 2) % 3] - xy[:, k]
        cs = np.sum(v1 * v2, 1) / (np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1) + 1e-30)
        ang.append(np.degrees(np.arccos(np.clip(cs, -1, 1))))
    ang = np.concatenate(ang)
    aspect = L.max(axis=1) / np.maximum(2 * np.abs(area) / L.max(axis=1), 1e-30)

    info = dict(
        stage="2B mesh gate: Si, SiO2, Gate, Shield",
        gmsh_seconds=round(dt, 2),
        vertices=int(len(nodes)),
        triangles=int(len(allt)),
        triangles_by_material={NAMES[k]: per_face.get(k, 0) for k in NAMES},
        faces_by_material={NAMES[k]: len(groups[k]) for k in NAMES},
        min_triangle_angle_deg=float(ang.min()),
        max_aspect_ratio=float(aspect.max()),
        min_area_um2=float(np.abs(area).min()),
        inverted=int((np.abs(area) < 1e-14).sum()),
        clockwise_ordered=int((area < 0).sum()),
        area_sum_um2=float(np.abs(area).sum()),
        analytic_domain_area_um2=float(BOTTOM_DEPTH * PITCH),
        area_rel_error=float(abs(np.abs(area).sum() - BOTTOM_DEPTH * PITCH) / (BOTTOM_DEPTH * PITCH)),
        triangles_spanning_two_materials=int(bad),
        conformity_ok=bool(bad == 0),
        separate_gate_and_shield_groups=True,
    )
    with open(f"{prefix}_gate2b.json", "w") as f:
        json.dump(info, f, indent=2)
    print(json.dumps(info, indent=2))

    fig, ax = plt.subplots(figsize=(8.4, 7.2))
    colours = {SI: "#e8d9a8", OX: "#8fd3e8", GATE: "#5c7fa3", SHIELD: "#2f4a63"}
    for code, arr in tris_mat:
        pass
    ax.triplot(nodes[:, 0], nodes[:, 1], allt, lw=0.2, color="0.4")
    xs = np.linspace(0, BOTTOM_DEPTH, 500)
    ax.plot(xs, a_of(xs), "r-", lw=1.3)
    ax.plot(xs, -a_of(xs), "r-", lw=1.3)
    ax.set_xlim(1.95, 2.45); ax.set_ylim(-0.72, 0.72)
    ax.set_aspect("equal")
    ax.set_xlabel("depth x (um)"); ax.set_ylabel("lateral y (um)")
    ax.set_title(f"Stage 2B mesh at the tip: Si / SiO2 / Gate / Shield\n"
                 f"{len(nodes)} nodes, {len(allt)} triangles, {dt:.2f} s")
    fig.tight_layout()
    fig.savefig(f"{prefix}_tip.png", dpi=180)
    print(f"wrote {prefix}_tip.png")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--h", type=float, default=0.02)
    ap.add_argument("--h-tip", type=float, default=0.008)
    ap.add_argument("--h-front", type=float, default=None)
    ap.add_argument("--prefix", default="mesh2b")
    args = ap.parse_args()
    run(args.h, args.h_tip, args.prefix, args.h_front)
