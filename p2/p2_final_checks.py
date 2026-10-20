#!/usr/bin/env python3
"""Final two P2.1 re-validation checks on the corrected solutions.

A. Classified interface audit on p2_r3_fixed.npz / p2_r4_fixed.npz:
   wall | arc proper (2.04<x<=2.14) | cap+corner (x>2.14)
   electrode-adjacent (all Si-element nodes pinned) reported separately
   near-zero pairs (D_scale < 1e-3 D_ref) excluded from the relative
   distribution only; every sample stays in the absolute integral.

B. Same-mesh P1/P2 probe comparison at the 84 archived probes, split into
   near-drain (x>2.38) and other, plus an explicit corner set around the
   79.51-degree corner at (2.15, +/-0.525).
"""
from __future__ import annotations
import json
import numpy as np
import p2_interface as PI
import p2_audit as A

EPS_SI, EPS_OX = A.EPS_SI, A.EPS_OX
DREF = 282.9836774371259          # frozen from p2_nearzero.json


# ---------------- A ----------------
def samples(fn):
    z = np.load(fn)
    nodes, tris, mat, phi, fixed = z["nodes"], z["tris"], z["mat"], z["phi"], z["fixed"]
    _, _, _, g1, g2, g3 = A.element_geom(nodes, tris)
    rows = []
    for ed in PI.build_edge_map(nodes, tris, mat):
        a, b, m = ed["nodes"]
        p0, pm, p1 = nodes[a], nodes[m], nodes[b]
        si_el, ox_el = ed["si"], ed["ox"]
        npin = int(fixed[tris[si_el]].sum())
        allpin = bool(np.all(fixed[tris[si_el]]))
        for t, w in zip(PI.G1, PI.GW):
            pt, dpt = PI.edge_points(p0, pm, p1, t)
            ds = np.linalg.norm(dpt)*w
            tn = dpt/np.linalg.norm(dpt)
            nrm = np.array([tn[1], -tn[0]])
            cen = nodes[tris[si_el][:3]].mean(axis=0)
            if np.dot(cen - pt, nrm) > 0:
                nrm = -nrm
            vals = []
            for e in (si_el, ox_el):
                P = nodes[tris[e]][:3]
                Mx = np.array([[P[1, 0]-P[0, 0], P[2, 0]-P[0, 0]],
                               [P[1, 1]-P[0, 1], P[2, 1]-P[0, 1]]])
                st = np.linalg.solve(Mx, pt - P[0])
                G = A.shape_grads(st[0], st[1], g1[e], g2[e], g3[e])
                vals.append(np.array([G[:, 0] @ phi[tris[e]], G[:, 1] @ phi[tris[e]]]))
            reg = "wall" if pt[0] <= 2.04 else ("arc" if pt[0] <= 2.14 else "cap")
            rows.append(dict(x=float(pt[0]), y=float(pt[1]), ds=float(ds),
                             d_si=EPS_SI*float(vals[0] @ nrm),
                             d_ox=EPS_OX*float(vals[1] @ nrm),
                             region=reg, npin=npin, allpin=allpin))
    return rows


def stats(rows):
    if not rows:
        return None
    w = np.array([r["ds"] for r in rows])
    ds = np.array([r["d_si"] for r in rows])
    do = np.array([r["d_ox"] for r in rows])
    den = np.sum(w*(np.abs(ds)+np.abs(do))/2)
    dscale = (np.abs(ds)+np.abs(do))/2
    nz = dscale < 1e-3*DREF
    rel = np.abs(ds-do)/np.maximum(2*dscale, 1e-30)
    keep = ~nz
    out = dict(n=len(rows),
               R_abs=float(np.sum(w*np.abs(ds-do))/den),
               R_signed=float(abs(np.sum(w*(ds-do)))/den),
               near_zero=int(nz.sum()),
               rel_median=float(np.median(rel[keep])) if keep.any() else None,
               rel_p95=float(np.percentile(rel[keep], 95)) if keep.any() else None,
               rel_max=float(np.max(rel[keep])) if keep.any() else None)
    if nz.any():
        xs = [r["x"] for r, m in zip(rows, nz) if m]
        ys = [r["y"] for r, m in zip(rows, nz) if m]
        out["near_zero_x"] = [float(min(xs)), float(max(xs))]
        out["near_zero_y"] = [float(min(ys)), float(max(ys))]
    return out


def part_a():
    out = {}
    for tag, fn in (("r3", "p2_r3_fixed.npz"), ("r4", "p2_r4_fixed.npz")):
        rows = samples(fn)
        out[tag] = {}
        print(f"\n[{tag}] classified interface audit  samples={len(rows)}")
        for reg in ("wall", "arc", "cap"):
            rr = [r for r in rows if r["region"] == reg]
            iface = [r for r in rr if not r["allpin"]]
            elec = [r for r in rr if r["allpin"]]
            s_all = stats(rr)
            s_if = stats(iface)
            out[tag][reg] = dict(all=s_all, interface=s_if,
                                 electrode_n=len(elec))
            if s_all:
                print(f"  {reg:5s} n={s_all['n']:4d} R_abs={s_all['R_abs']:.6f} "
                      f"R_sgn={s_all['R_signed']:.6f} nz={s_all['near_zero']:3d} | "
                      f"iface n={s_if['n']:4d} R_abs={s_if['R_abs']:.6f} "
                      f"R_sgn={s_if['R_signed']:.6f} med={s_if['rel_median']:.5f} "
                      f"p95={s_if['rel_p95']:.5f} max={s_if['rel_max']:.4f} "
                      f"| elec={len(elec)}")
        s_all = stats(rows); s_if = stats([r for r in rows if not r["allpin"]])
        out[tag]["all"] = dict(all=s_all, interface=s_if)
        print(f"  all   n={s_all['n']:4d} R_abs={s_all['R_abs']:.6f} "
              f"R_sgn={s_all['R_signed']:.6f} | iface n={s_if['n']:4d} "
              f"R_abs={s_if['R_abs']:.6f} R_sgn={s_if['R_signed']:.6f}")
    json.dump(out, open("p2_classified_fixed.json", "w"), indent=2)


# ---------------- B ----------------
def locate(nodes, tris, pts):
    Aa = nodes[tris[:, 0]]; Bb = nodes[tris[:, 1]]; Cc = nodes[tris[:, 2]]
    det = ((Bb[:, 0]-Aa[:, 0])*(Cc[:, 1]-Aa[:, 1]) -
           (Cc[:, 0]-Aa[:, 0])*(Bb[:, 1]-Aa[:, 1]))
    el = np.full(len(pts), -1, int); bar = np.zeros((len(pts), 2))
    for k, p in enumerate(pts):
        L1 = ((Bb[:, 1]-Cc[:, 1])*(p[0]-Cc[:, 0]) +
              (Cc[:, 0]-Bb[:, 0])*(p[1]-Cc[:, 1]))/det
        L2 = ((Cc[:, 1]-Aa[:, 1])*(p[0]-Cc[:, 0]) +
              (Aa[:, 0]-Cc[:, 0])*(p[1]-Cc[:, 1]))/det
        L3 = 1-L1-L2
        m = (L1 >= -1e-9) & (L2 >= -1e-9) & (L3 >= -1e-9)
        if m.any():
            j = int(np.flatnonzero(m)[0])
            el[k] = j; bar[k] = (L1[j], L2[j])
    return el, bar


def eval_p1(nodes, tris, phi, pts):
    el, bar = locate(nodes, tris, pts)
    out = np.full(len(pts), np.nan)
    for k in range(len(pts)):
        if el[k] < 0:
            continue
        L1, L2 = bar[k]
        ids = tris[el[k]]
        out[k] = L1*phi[ids[0]] + L2*phi[ids[1]] + (1-L1-L2)*phi[ids[2]]
    return out


def eval_p2(nodes, tris, phi, pts):
    el, bar = locate(nodes, tris, pts)
    out = np.full(len(pts), np.nan)
    for k in range(len(pts)):
        if el[k] < 0:
            continue
        L1, L2 = bar[k]
        out[k] = float(np.dot(phi[tris[el[k]]], A.shape_vals(L1, L2)))
    return out


def part_b():
    p1 = np.load("p1_on_p2mesh.npz")
    p2 = np.load("p2_r3_variant.npz")
    xs = np.linspace(2.25, 2.415, 12)
    ys = np.array([-0.30, -0.20, -0.10, 0.0, 0.10, 0.20, 0.30])
    pts84 = np.array([[x, y] for x in xs for y in ys])

    corner = []
    for cy in (0.525, -0.525):
        for r in (0.005, 0.02, 0.04):
            for th in np.linspace(0, 2*np.pi, 12, endpoint=False):
                corner.append([2.15 + r*np.cos(th), cy + r*np.sin(th)])
    corner = np.array(corner)

    f1 = eval_p1(p1["nodes"], p1["tris"], p1["phi"], pts84)
    f2 = eval_p2(p2["nodes"], p2["tris"], p2["phi"], pts84)
    ok = ~np.isnan(f1) & ~np.isnan(f2)
    d = f2[ok] - f1[ok]
    xo = np.repeat(xs, len(ys))[ok]
    groups = {"near_drain_x>2.38": xo > 2.38, "archived_other": xo <= 2.38}
    print("\n[B] same-mesh probes (corrected P2)")
    out = {}
    for name, m in groups.items():
        out[name] = dict(n=int(m.sum()),
                         rms=float(np.sqrt(np.mean(d[m]**2))) if m.any() else None,
                         max_abs=float(np.abs(d[m]).max()) if m.any() else None,
                         mean=float(d[m].mean()) if m.any() else None)
        print(f"  {name:18s} n={out[name]['n']:3d} rms={out[name]['rms']:.6f} "
              f"max={out[name]['max_abs']:.6f} mean={out[name]['mean']:+.6f}")

    g1 = eval_p1(p1["nodes"], p1["tris"], p1["phi"], corner)
    g2 = eval_p2(p2["nodes"], p2["tris"], p2["phi"], corner)
    okc = ~np.isnan(g1) & ~np.isnan(g2)
    dc = g2[okc] - g1[okc]
    out["corner"] = dict(n=int(okc.sum()),
                         rms=float(np.sqrt(np.mean(dc**2))) if okc.any() else None,
                         max_abs=float(np.abs(dc).max()) if okc.any() else None,
                         mean=float(dc.mean()) if okc.any() else None)
    print(f"  {'corner (2.15,.525)':18s} n={out['corner']['n']:3d} "
          f"rms={out['corner']['rms']:.6f} max={out['corner']['max_abs']:.6f} "
          f"mean={out['corner']['mean']:+.6f}")
    out["all84"] = dict(n=int(ok.sum()), rms=float(np.sqrt(np.mean(d**2))),
                        max_abs=float(np.abs(d).max()), mean=float(d.mean()))
    print(f"  {'all archived':18s} n={out['all84']['n']:3d} rms={out['all84']['rms']:.6f} "
          f"max={out['all84']['max_abs']:.6f} mean={out['all84']['mean']:+.6f}")
    json.dump(out, open("p2_probe_fixed.json", "w"), indent=2)


if __name__ == "__main__":
    part_a()
    part_b()
    print("\nwrote p2_classified_fixed.json and p2_probe_fixed.json")
