#!/usr/bin/env python3
"""Final P2.1 interface table on both resolutions, with the electrode-adjacent
class separated from the no-free-charge class.

Classification (frozen, applied identically to both grids):
  electrode_adjacent : the adjacent Si element has all six nodes Dirichlet-pinned
                       (equipotential contact) -> eps_Si grad(phi).n == 0 there
                       by construction, and the no-free-charge condition does
                       not apply; flux legitimately terminates on the electrode.
  near_zero          : D_scale < 1e-3 * D_ref  (D_ref frozen at 282.983677)
  interface          : everything else
All classes remain in the integrated R_abs; only the relative distribution
excludes near_zero, which is reported separately with its count.
"""
from __future__ import annotations
import json
import numpy as np
import p2_interface as PI
import p2_audit as A

EPS_SI, EPS_OX = A.EPS_SI, A.EPS_OX


def samples(fn, dref):
    z = np.load(fn)
    nodes, tris, mat, phi, fixed = z["nodes"], z["tris"], z["mat"], z["phi"], z["fixed"]
    p, det, area, g1, g2, g3 = A.element_geom(nodes, tris)
    rows = []
    for ed in PI.build_edge_map(nodes, tris, mat):
        a, b, m = ed["nodes"]
        p0, pm, p1 = nodes[a], nodes[m], nodes[b]
        pinned = bool(np.all(fixed[tris[ed["si"]]]))
        for t, w in zip(PI.G1, PI.GW):
            pt, dpt = PI.edge_points(p0, pm, p1, t)
            ds = np.linalg.norm(dpt) * w
            tn = dpt / np.linalg.norm(dpt)
            nrm = np.array([tn[1], -tn[0]])
            cen = nodes[tris[ed["si"]][:3]].mean(axis=0)
            if np.dot(cen - pt, nrm) > 0:
                nrm = -nrm
            vals = []
            for e in (ed["si"], ed["ox"]):
                P = nodes[tris[e]][:3]
                Mx = np.array([[P[1,0]-P[0,0], P[2,0]-P[0,0]],
                               [P[1,1]-P[0,1], P[2,1]-P[0,1]]])
                st = np.linalg.solve(Mx, pt - P[0])
                vals.append(PI.grad_at(nodes, tris, phi, e, st[0], st[1], g1, g2, g3))
            rows.append(dict(x=float(pt[0]), y=float(pt[1]), ds=float(ds),
                             d_si=EPS_SI*float(vals[0] @ nrm),
                             d_ox=EPS_OX*float(vals[1] @ nrm),
                             region="arc" if pt[0] > 2.04 else "wall",
                             pinned=pinned))
    return rows


def summarize(rows, label, dref):
    if not rows:
        return dict(label=label, n=0)
    w = np.array([r["ds"] for r in rows])
    ds = np.array([r["d_si"] for r in rows])
    do = np.array([r["d_ox"] for r in rows])
    dscale = (np.abs(ds)+np.abs(do))/2
    den = np.sum(w*(np.abs(ds)+np.abs(do))/2)
    rel = np.abs(ds-do)/np.maximum(2*dscale, 1e-30)
    nz = dscale < 1e-3*dref
    keep = ~nz
    return dict(label=label, n=len(rows),
                R_abs=float(np.sum(w*np.abs(ds-do))/den),
                R_signed=float(abs(np.sum(w*(ds-do)))/den),
                rel_median=float(np.median(rel[keep])) if keep.any() else None,
                rel_mean=float(np.mean(rel[keep])) if keep.any() else None,
                rel_p95=float(np.percentile(rel[keep], 95)) if keep.any() else None,
                rel_max=float(np.max(rel[keep])) if keep.any() else None,
                near_zero=int(nz.sum()),
                opp_sign=int((np.sign(ds) != np.sign(do)).sum()),
                mean_d_si=float(ds.mean()), mean_d_ox=float(do.mean()))


def main():
    dref = json.load(open("p2_nearzero.json"))["reference"]["d_ref"]
    out = {"d_ref": dref, "resolutions": {}}
    for tag, fn, h in (("r3", "p2_r3_solution.npz", 0.014),
                       ("r4", "p2_r4_solution.npz", 0.010)):
        rows = samples(fn, dref)
        wall = [r for r in rows]
        res = dict(h=h, total=len(rows))
        for reg in ("wall", "arc"):
            rr = [r for r in rows if r["region"] == reg]
            res[reg] = dict(
                all=summarize(rr, f"{tag}/{reg}/all", dref),
                interface=summarize([r for r in rr if not r["pinned"]], f"{tag}/{reg}/interface", dref),
                electrode=summarize([r for r in rr if r["pinned"]], f"{tag}/{reg}/electrode", dref),
            )
        out["resolutions"][tag] = res
        for reg in ("wall", "arc"):
            for cls in ("all", "interface", "electrode"):
                s = res[reg][cls]
                if s.get("n"):
                    print(f"{tag} {reg:4s} {cls:9s} n={s['n']:5d} "
                          f"R_abs={s['R_abs']:.6f} R_sgn={s['R_signed']:.6f} "
                          f"med={s['rel_median']:.5f} p95={s['rel_p95']:.5f} "
                          f"max={s['rel_max']:.4f} nz={s['near_zero']:3d} "
                          f"Dsi={s['mean_d_si']:8.2f} Dox={s['mean_d_ox']:8.2f}")
    with open("p2_final_table.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print("wrote p2_final_table.json")


if __name__ == "__main__":
    main()
