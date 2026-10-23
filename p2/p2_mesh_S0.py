#!/usr/bin/env python3
"""S0 pilot mesh generator (L0 / L1), global scaling.

h_k(x,y) = r^k * h0(x,y),  r = 2^-1/2
  band core        0.0040 * r^k
  background       0.0085 * r^k
  tip              0.0034 * r^k
  band distances   constant (0.04 / 0.06 um)  -- widths do not scale

PILOT ONLY: generates a mesh and measures it. No assembly, no solve.
"""
from __future__ import annotations
import json
import sys
import numpy as np
import gmsh
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import mesh_stage2b as M
import p2_audit as A

R = 2.0**-0.5
H0 = dict(band=0.0040, background=0.0085, tip=0.0034)
DIST_MIN, DIST_MAX = 0.04, 0.06
N_ARC = 65          # canonical chord count per circular arc (64 chords)


def build(level, outdir):
    rk = R**level
    hb = H0["band"]*rk
    hg = H0["background"]*rk
    ht = H0["tip"]*rk
    prefix = f"S0_L{level}"
    print(f"[{prefix}] r^k={rk:.6f} band={hb:.6f} background={hg:.6f} tip={ht:.6f}",
          flush=True)

    band = np.load("p2_band_contour.npz")["pts"]
    sel = band[::4]

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 1)
    groups = M.build(hg, ht)

    occ = gmsh.model.occ
    tags = [occ.addPoint(float(x), float(y), 0.0, hb) for x, y in sel]
    spl = occ.addSpline(tags)
    occ.synchronize()

    # Canonical geometry: force a fixed chord count on every circular arc so the
    # polygonal domain is identical at all levels (CP-M1a). Query by curve type
    # because OCC fragmentation renumbers entities.
    arcs = [abs(t) for d, t in gmsh.model.getEntities(1)
            if gmsh.model.getType(1, abs(t)) == "Circle"]
    for t in arcs:
        gmsh.model.mesh.setTransfiniteCurve(t, N_ARC)
    print(f"[{prefix}] canonical arcs: {len(arcs)} x {N_ARC} nodes", flush=True)

    fd = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(fd, "CurvesList", [spl])
    gmsh.model.mesh.field.setNumber(fd, "Sampling", 200)
    ft = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(ft, "InField", fd)
    gmsh.model.mesh.field.setNumber(ft, "SizeMin", hb)
    gmsh.model.mesh.field.setNumber(ft, "SizeMax", hg)
    gmsh.model.mesh.field.setNumber(ft, "DistMin", DIST_MIN)
    gmsh.model.mesh.field.setNumber(ft, "DistMax", DIST_MAX)
    gmsh.model.mesh.field.setAsBackgroundMesh(ft)

    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.SecondOrderLinear", 1)
    gmsh.option.setNumber("Mesh.MeshSizeMin", min(ht, hb))
    gmsh.option.setNumber("Mesh.MeshSizeMax", hg)
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
    np.savez_compressed(f"{outdir}/{prefix}_p2mesh.npz", nodes=nodes, tris=tris, mat=mat)

    P = nodes[tris][:, :3]
    det = ((P[:, 1, 0]-P[:, 0, 0])*(P[:, 2, 1]-P[:, 0, 1]) -
           (P[:, 2, 0]-P[:, 0, 0])*(P[:, 1, 1]-P[:, 0, 1]))
    area = 0.5*np.abs(det)
    si = mat == M.SI
    cx = P[:, :, 0].mean(1); cy = P[:, :, 1].mean(1)
    charged = si & (cx > A.BODY_BOTTOM + 1e-12)
    from p2_charge_mask_audit import clip_area
    xmin = P[:, :, 0].min(1); xmax = P[:, :, 0].max(1)
    exact = area[si & (xmin > A.BODY_BOTTOM)].sum()
    for k in np.flatnonzero(si & (xmin <= A.BODY_BOTTOM) & (xmax > A.BODY_BOTTOM)):
        exact += clip_area([tuple(q) for q in P[k]], A.BODY_BOTTOM)

    h = np.sqrt(4*area/np.sqrt(3))
    from scipy.spatial import cKDTree
    d, _ = cKDTree(band).query(np.column_stack([cx, cy]))
    regs = {"band<=0.04": d <= 0.04,
            "trans0.04-0.06": (d > 0.04) & (d <= 0.06),
            "tip": (cx > 2.30) & (cx < 2.42) & (np.abs(cy) < 0.02),
            "bodycut": np.abs(cx - A.BODY_BOTTOM) <= 0.05,
            "background": d > 0.15}
    stats = {}
    for k, m in regs.items():
        m = m & si
        if m.sum() < 5:
            continue
        stats[k] = dict(n=int(m.sum()), p50=float(np.median(h[m])),
                        p90=float(np.percentile(h[m], 90)),
                        p95=float(np.percentile(h[m], 95)))
    edges = np.stack([np.linalg.norm(P[:, 1]-P[:, 0], axis=1),
                      np.linalg.norm(P[:, 2]-P[:, 1], axis=1),
                      np.linalg.norm(P[:, 0]-P[:, 2], axis=1)], 1)
    emax = edges.max(1); hmin = 2*area/np.maximum(emax, 1e-30)
    rec = dict(level=level, rk=rk, n_arc=N_ARC, n_arcs=len(arcs), band=hb, background=hg, tip=ht,
               nodes=int(len(nodes)), elements=int(len(tris)), counts=counts,
               si_area=float(area[si].sum()), exact_charged=float(exact),
               centroid_charged=float(area[charged].sum()),
               aspect_max=float((emax/np.maximum(hmin, 1e-30)).max()),
               regions=stats)
    json.dump(rec, open(f"{outdir}/{prefix}_mesh.json", "w"), indent=2)
    print(f"[{prefix}] nodes={len(nodes)} elems={len(tris)} "
          f"Si={rec['si_area']:.12f} exact={rec['exact_charged']:.12f}", flush=True)
    return rec


if __name__ == "__main__":
    out = "."
    out = os.path.join("..","p2_3_spatial_S0","levels")
    os.makedirs(out, exist_ok=True)
    recs = [build(0, out), build(1, out)]
    json.dump(recs, open(f"{out}/pilot_summary.json", "w"), indent=2)
    print("wrote pilot_summary.json")
