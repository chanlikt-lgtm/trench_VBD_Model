#!/usr/bin/env python3
"""Stage 2A: geometry + mesh validation only. No electrostatics.

Geometry is built from explicit OCC curves that are SHARED between the silicon
and oxide faces (no embedded splines, no boundary-to-boundary spline):

  mesa profile (upper, x = depth, y = lateral):
      (0, 0.340) --line--> (2.05, 0.425) --arc r=0.10 about (2.05,0.525)-->
      (2.15, 0.525) --line--> (2.42, 0.575)
  mirrored for the lower half; the two chains are closed at x=0 and x=BOTTOM.

The mesa half width is capped at A_MAX = P/2 - t_ox = 0.575 um, so the oxide
keeps its full 0.10 um everywhere including the tip. The oxide face therefore
closes to a corner at (2.42, +/-0.675) rather than collapsing to zero width.
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
X_ARC_END = TIP_START + R_TIP          # 2.15
A_ARC_END = A_BOT + R_TIP              # 0.525
BOTTOM_DEPTH = 2.42
A_MAX = PITCH / 2 - 0.10               # 0.575
OXIDE_T = 0.10
Y_EDGE = PITCH / 2


def build_geometry(h_global, h_tip):
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)
    gmsh.model.add("trench2a")
    occ = gmsh.model.occ

    rect = occ.addRectangle(0.0, -Y_EDGE, 0.0, BOTTOM_DEPTH, PITCH)

    # --- profile points, upper half (shared between Si and oxide) ---
    p0 = occ.addPoint(0.0, A_TOP, 0.0, h_global)
    p1 = occ.addPoint(TIP_START, A_BOT, 0.0, h_tip)
    pc = occ.addPoint(TIP_START, A_ARC_END, 0.0, h_tip)     # arc centre
    p2 = occ.addPoint(X_ARC_END, A_ARC_END, 0.0, h_tip)
    p3 = occ.addPoint(BOTTOM_DEPTH, A_MAX, 0.0, h_global)

    up_l1 = occ.addLine(p0, p1)
    up_arc = occ.addCircleArc(p1, pc, p2)
    up_l2 = occ.addLine(p2, p3)

    # mirrored lower half
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

    loop = occ.addCurveLoop([up_l1, up_arc, up_l2, right, lo_l2, lo_arc, lo_l1, left])
    si = occ.addPlaneSurface([loop])
    occ.synchronize()

    # fragment: rectangle minus the mesa, sharing the profile curves
    out, _ = occ.fragment([(2, rect)], [(2, si)])
    occ.synchronize()

    faces = [t for d, t in out if d == 2]
    com = {t: np.array(gmsh.model.occ.getCenterOfMass(2, t)) for t in faces}

    # classify by centroid: inside the mesa -> Si, else oxide
    si_faces, ox_faces = [], []
    for t in faces:
        cx, cy = com[t][0], com[t][1]
        a = half_width_scalar(cx)
        (si_faces if abs(cy) < a else ox_faces).append(t)

    # physical groups: shared interface curve = the profile curves (one side)
    gmsh.model.addPhysicalGroup(2, si_faces, name="Si")
    gmsh.model.addPhysicalGroup(2, ox_faces, name="SiO2")
    gmsh.model.addPhysicalGroup(1, [up_l1, up_arc, up_l2], name="interface_upper")
    gmsh.model.addPhysicalGroup(1, [lo_l1, lo_arc, lo_l2], name="interface_lower")

    print(f"[geom] {len(faces)} faces: {len(si_faces)} Si, {len(ox_faces)} oxide", flush=True)
    return dict(si=si_faces, ox=ox_faces, rect=rect)


def half_width_scalar(x):
    if x <= TIP_START:
        return A_TOP + (A_BOT - A_TOP) * x / TIP_START
    if x <= X_ARC_END:
        u = (x - TIP_START) / R_TIP
        return A_BOT + R_TIP - R_TIP * np.sqrt(max(0.0, 1.0 - u * u))
    t = (x - X_ARC_END) / (BOTTOM_DEPTH - X_ARC_END)
    return A_ARC_END + (A_MAX - A_ARC_END) * t


def mesh_and_report(h_global, h_tip, prefix):
    t0 = time.time()
    tags = build_geometry(h_global, h_tip)

    gmsh.option.setNumber("Mesh.MeshSizeMin", h_tip)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h_global)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.model.mesh.generate(2)
    dt = time.time() - t0
    print(f"[mesh] generated in {dt:.2f} s", flush=True)

    ntags, ncoords, _ = gmsh.model.mesh.getNodes()
    nodes = np.column_stack([ncoords[0::3], ncoords[1::3]])
    tag2i = {int(t): i for i, t in enumerate(ntags)}

    per_face = {}
    tris_all = []
    for name, flist in (("Si", tags["si"]), ("SiO2", tags["ox"])):
        cnt = 0
        for f in flist:
            etypes, etags, enodes = gmsh.model.mesh.getElements(2, f)
            for et, en in zip(etypes, enodes):
                if et != 2:
                    continue
                arr = np.asarray(en, dtype=int).reshape(-1, 3)
                tris_all.append(np.vectorize(tag2i.get)(arr))
                cnt += len(arr)
        per_face[name] = cnt

    tris = np.vstack(tris_all) if tris_all else np.zeros((0, 3), int)
    gmsh.finalize()

    xy = nodes[tris]
    e1 = xy[:, 1] - xy[:, 0]; e2 = xy[:, 2] - xy[:, 1]; e3 = xy[:, 0] - xy[:, 2]
    L = np.stack([np.linalg.norm(e1, axis=1), np.linalg.norm(e2, axis=1),
                  np.linalg.norm(e3, axis=1)], 1)
    cr = e1[:, 0] * (xy[:, 2, 1] - xy[:, 0, 1]) - e1[:, 1] * (xy[:, 2, 0] - xy[:, 0, 0])
    area = 0.5 * np.abs(cr)
    ang = []
    for k in range(3):
        v1 = xy[:, (k + 1) % 3] - xy[:, k]
        v2 = xy[:, (k + 2) % 3] - xy[:, k]
        cs = np.sum(v1 * v2, 1) / (np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1) + 1e-30)
        ang.append(np.degrees(np.arccos(np.clip(cs, -1, 1))))
    ang = np.concatenate(ang)

    min_area = area.min()
    min_alt = 2 * area / L.max(axis=1)
    aspect = L.max(axis=1) / np.maximum(min_alt, 1e-30)
    inverted = int((area <= 0).sum())

    tip = area > 0
    near_tip = (xy[:, :, 0].mean(axis=1) > TIP_START - 0.12) & \
               (xy[:, :, 0].mean(axis=1) < BOTTOM_DEPTH + 0.01) & \
               (np.abs(xy[:, :, 1].mean(axis=1)) > 0.35)

    info = dict(
        stage="2A geometry + mesh only",
        gmsh_seconds=round(dt, 2),
        vertices=int(len(nodes)), triangles=int(len(tris)),
        triangles_Si=per_face.get("Si", 0), triangles_SiO2=per_face.get("SiO2", 0),
        min_triangle_angle_deg=float(ang.min()),
        mean_triangle_angle_deg=float(ang.mean()),
        max_aspect_ratio=float(aspect.max()),
        max_aspect_ratio_near_tip=float(aspect[near_tip].max()) if near_tip.any() else None,
        min_area_um2=float(min_area),
        inverted_or_zero_area=int(inverted),
        tip_radius_um=R_TIP, a_max_um=A_MAX, oxide_um=OXIDE_T,
        oxide_thickness_at_tip_um=float(Y_EDGE - A_MAX),
    )
    with open(f"{prefix}_mesh_gate.json", "w") as f:
        json.dump(info, f, indent=2)
    print(json.dumps(info, indent=2))

    # tip close-up
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.triplot(nodes[:, 0], nodes[:, 1], tris, lw=0.25, color="0.35")
    xs = np.linspace(0, BOTTOM_DEPTH, 600)
    ax.plot(xs, [half_width_scalar(v) for v in xs], "r-", lw=1.4, label="Si/SiO$_2$ interface")
    ax.plot(xs, [-half_width_scalar(v) for v in xs], "r-", lw=1.4)
    ax.set_xlim(1.95, 2.45); ax.set_ylim(-0.72, 0.72)
    ax.set_aspect("equal")
    ax.set_xlabel("depth x (um)"); ax.set_ylabel("lateral y (um)")
    ax.set_title(f"Stage 2A mesh at the tip (r_tip = {R_TIP} um)\n"
                 f"{len(nodes)} nodes, {len(tris)} triangles, {dt:.2f} s")
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{prefix}_tip.png", dpi=180)
    print(f"wrote {prefix}_tip.png")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--h", type=float, default=0.02)
    ap.add_argument("--h-tip", type=float, default=0.008)
    ap.add_argument("--prefix", default="mesh2a")
    a = ap.parse_args()
    mesh_and_report(a.h, a.h_tip, a.prefix)
