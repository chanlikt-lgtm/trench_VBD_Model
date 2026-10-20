#!/usr/bin/env python3
"""L1 mesh: intermediate background with graded local refinement.

Windows (physics unchanged):
  A: 2.32 <= x <= 2.37, |y| <= 0.13
  B: 2.395 <= x <= 2.420, |y| <= 0.24
"""
from __future__ import annotations
import json
import numpy as np
import gmsh
import sys, os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import mesh_stage2b as M

H, H_TIP, H_LOCAL = 0.0085, 0.0034, 0.0040
BOXES = [
    (2.32, 2.37, -0.13, 0.13, H_LOCAL),
    (2.395, 2.420, -0.24, 0.24, H_LOCAL),
]
THICK = 0.012
PREFIX = "p2_L1"


def main():
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 1)
    groups = M.build(H, H_TIP)
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.SecondOrderLinear", 1)
    gmsh.option.setNumber("Mesh.MeshSizeMin", min(H_TIP, H_LOCAL))
    gmsh.option.setNumber("Mesh.MeshSizeMax", H)

    f = gmsh.model.mesh.field
    ids = []
    for (x0, x1, y0, y1, hv) in BOXES:
        b = f.add("Box")
        f.setNumber(b, "VIn", hv)
        f.setNumber(b, "VOut", H)
        f.setNumber(b, "XMin", x0); f.setNumber(b, "XMax", x1)
        f.setNumber(b, "YMin", y0); f.setNumber(b, "YMax", y1)
        f.setNumber(b, "Thickness", THICK)
        ids.append(b)
    mn = f.add("Min")
    f.setNumbers(mn, "FieldsList", ids)
    f.setAsBackgroundMesh(mn)

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
    print(f"[L1] nodes={len(nodes)} elements={len(tris)} counts={counts}")

    # quality
    p = nodes[tris][:, :3]
    det = ((p[:, 1, 0]-p[:, 0, 0])*(p[:, 2, 1]-p[:, 0, 1]) -
           (p[:, 2, 0]-p[:, 0, 0])*(p[:, 1, 1]-p[:, 0, 1]))
    area = 0.5*np.abs(det)
    edges = np.stack([np.linalg.norm(p[:, 1]-p[:, 0], axis=1),
                      np.linalg.norm(p[:, 2]-p[:, 1], axis=1),
                      np.linalg.norm(p[:, 0]-p[:, 2], axis=1)], 1)
    emax = edges.max(axis=1)
    hmin = 2*area/np.maximum(emax, 1e-30)
    aspect = emax/np.maximum(hmin, 1e-30)
    ang = []
    for k in range(3):
        v1 = p[:, (k+1) % 3]-p[:, k]; v2 = p[:, (k+2) % 3]-p[:, k]
        cs = np.sum(v1*v2, 1)/(np.linalg.norm(v1, axis=1)*np.linalg.norm(v2, axis=1)+1e-30)
        ang.append(np.degrees(np.arccos(np.clip(cs, -1, 1))))
    ang = np.concatenate(ang)
    cxs = p[:, :, 0].mean(axis=1)
    si = mat == M.SI
    charged = si & (cxs > 0.54)
    print(f"  min area={area.min():.3e} max aspect={aspect.max():.2f} "
          f"min angle={ang.min():.2f} deg")
    print(f"  Si area={area[si].sum():.9f} charged area={area[charged].sum():.9f}")
    json.dump(dict(nodes=int(len(nodes)), elements=int(len(tris)),
                   counts={str(k): v for k, v in counts.items()},
                   min_area=float(area.min()), max_aspect=float(aspect.max()),
                   min_angle=float(ang.min()),
                   si_area=float(area[si].sum()),
                   charged_area=float(area[charged].sum())),
              open(f"{PREFIX}_mesh.json", "w"), indent=2)


if __name__ == "__main__":
    main()
