#!/usr/bin/env python3
"""Level-1 controlled-band mesh.

Band geometry is frozen in p2_band_contour.npz (L1b, eps = 1e-8).
Size policy:
    distance <= 0.04 um  -> h_band = 0.0040
    distance >= 0.06 um  -> h_global = 0.0085   (linear threshold between)
Background, geometry, material tags, charge rule and BCs unchanged.
"""
from __future__ import annotations
import json
import numpy as np
import gmsh
import sys, os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import mesh_stage2b as M
import p2_audit as A

H, H_TIP, H_BAND = 0.0085, 0.0034, 0.0040
DIST_MIN, DIST_MAX = 0.04, 0.06
PREFIX = "p2_B1"
REG = {"A": (2.32, 2.37, -0.13, 0.13), "B": (2.395, 2.420, -0.24, 0.24)}


def main():
    band = np.load("p2_band_contour.npz")["pts"]
    sel = band[::4]                      # ~156 spline points
    print(f"[B1] band points {len(band)} -> spline {len(sel)}", flush=True)

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 1)
    groups = M.build(H, H_TIP)

    occ = gmsh.model.occ
    tags = [occ.addPoint(float(x), float(y), 0.0, H_BAND) for x, y in sel]
    spl = occ.addSpline(tags)
    occ.synchronize()

    fd = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(fd, "CurvesList", [spl])
    gmsh.model.mesh.field.setNumber(fd, "Sampling", 200)

    ft = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(ft, "InField", fd)
    gmsh.model.mesh.field.setNumber(ft, "SizeMin", H_BAND)
    gmsh.model.mesh.field.setNumber(ft, "SizeMax", H)
    gmsh.model.mesh.field.setNumber(ft, "DistMin", DIST_MIN)
    gmsh.model.mesh.field.setNumber(ft, "DistMax", DIST_MAX)
    gmsh.model.mesh.field.setAsBackgroundMesh(ft)

    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.SecondOrderLinear", 1)
    gmsh.option.setNumber("Mesh.MeshSizeMin", H_TIP)
    gmsh.option.setNumber("Mesh.MeshSizeMax", H)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    gmsh.model.mesh.generate(2)
    ntags, ncoords, _ = gmsh.model.mesh.getNodes()
    nodes = np.column_stack([ncoords[0::3], ncoords[1::3]])
    tag2i = {int(t): i for i, t in enumerate(ntags)}

    parts, mats, counts = {}, {}, {}
    for code, lst in groups.items():
        acc = []
        for face in lst:
            etypes, _, enodes = gmsh.model.mesh.getElements(2, face)
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
    np.savez_compressed(f"{PREFIX}_p2mesh.npz", nodes=nodes, tris=tris, mat=mat)
    print(f"[B1] nodes={len(nodes)} elements={len(tris)} counts={counts}", flush=True)

    p = nodes[tris][:, :3]
    det = ((p[:, 1, 0]-p[:, 0, 0])*(p[:, 2, 1]-p[:, 0, 1]) -
           (p[:, 2, 0]-p[:, 0, 0])*(p[:, 1, 1]-p[:, 0, 1]))
    area = 0.5*np.abs(det)
    edges = np.stack([np.linalg.norm(p[:, 1]-p[:, 0], axis=1),
                      np.linalg.norm(p[:, 2]-p[:, 1], axis=1),
                      np.linalg.norm(p[:, 0]-p[:, 2], axis=1)], 1)
    emax = edges.max(1); hmin = 2*area/np.maximum(emax, 1e-30)
    aspect = emax/np.maximum(hmin, 1e-30)
    ang = []
    for k in range(3):
        v1 = p[:, (k+1) % 3]-p[:, k]; v2 = p[:, (k+2) % 3]-p[:, k]
        cs = np.sum(v1*v2, 1)/(np.linalg.norm(v1, axis=1)*np.linalg.norm(v2, axis=1)+1e-30)
        ang.append(np.degrees(np.arccos(np.clip(cs, -1, 1))))
    ang = np.concatenate(ang)
    si = mat == M.SI
    cxs = p[:, :, 0].mean(1)
    charged = si & (cxs > A.BODY_BOTTOM + 1e-12)

    h = np.sqrt(4*area/np.sqrt(3))
    print(f"  min area={area.min():.3e} max aspect={aspect.max():.2f} "
          f"min angle={ang.min():.2f} deg")
    print(f"  Si area={area[si].sum():.12f} centroid charged={area[charged].sum():.12f}")

    from p2_charge_mask_audit import clip_area
    P = nodes[tris][:, :3]
    xmin = P[:, :, 0].min(1); xmax = P[:, :, 0].max(1)
    stretch = si & (xmin <= A.BODY_BOTTOM) & (xmax > A.BODY_BOTTOM)
    exact = area[si & (xmin > A.BODY_BOTTOM)].sum()
    for k in np.flatnonzero(stretch):
        exact += clip_area([tuple(q) for q in P[k]], A.BODY_BOTTOM)
    print(f"  exact clipped charged={exact:.12f}  straddle={int(stretch.sum())}")

    stats = {}
    for tag, rect in REG.items():
        m = (cxs >= rect[0]) & (cxs <= rect[1]) & (p[:, :, 1].mean(1) >= rect[2]) & \
            (p[:, :, 1].mean(1) <= rect[3])
        if m.sum() > 5:
            stats[tag] = dict(n=int(m.sum()), h_med=float(np.median(h[m])))
    strip = (cxs > 2.37) & (cxs < 2.395)
    if strip.sum() > 5:
        stats["strip"] = dict(n=int(strip.sum()), h_med=float(np.median(h[strip])))
    print("  window medians: " + "  ".join(f"{k}:{v['h_med']:.5f}" for k, v in stats.items()))

    json.dump(dict(nodes=int(len(nodes)), elements=int(len(tris)),
                   h_band=H_BAND, h_global=H, dist_min=DIST_MIN, dist_max=DIST_MAX,
                   counts={str(k): v for k, v in counts.items()},
                   min_area=float(area.min()), max_aspect=float(aspect.max()),
                   min_angle=float(ang.min()),
                   si_area=float(area[si].sum()),
                   centroid_charged=float(area[charged].sum()),
                   exact_charged=float(exact),
                   windows={k: v for k, v in stats.items()}),
              open(f"{PREFIX}_mesh.json", "w"), indent=2)
    print("wrote p2_B1_p2mesh.npz / p2_B1_mesh.json")


if __name__ == "__main__":
    main()
