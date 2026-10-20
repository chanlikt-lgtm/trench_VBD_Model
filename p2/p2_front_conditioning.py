#!/usr/bin/env python3
"""Front-conditioning diagnostic for the P2.3 neutral front.

For cuts x = 2.35 and 2.41, on the five F7-certified solutions:
  - main interval endpoints at eps = 1e-8, 1e-6, 1e-5
  - local u(y) and du/dy from the full P2 interpolant at each endpoint
  - conditioning estimate |dy| ~ |du| / |du/dy|
  - maximum u inside each gap between adjacent intervals
  - cross-mesh u differences on a common y grid

Node ordering is detected automatically (P2 mid-side permutation).
"""
from __future__ import annotations
import json
import itertools
import numpy as np

from p2_exact_metrics import cut_intervals
from p2_topology_clean import clean

MESHES = [
    ("r3",  "p2_p3_r3_polished.npz"),
    ("r4",  "p2_p3_r4_polished.npz"),
    ("mid", "p2_p3_rmid_polished.npz"),
    ("L1b", "p2_L1b_scaled.npz"),
    ("L1",  "p2_L1_scaled.npz"),
]
CUTS = (2.35, 2.41)
EPSS = (1e-8, 1e-6, 1e-5)
YG = np.linspace(-0.34, 0.34, 681)


def detect_order(nodes, tris):
    """Return the permutation mapping local mid-side slots 3,4,5 to edges."""
    p = nodes[tris]
    v = p[:, :3]
    mids = p[:, 3:]
    best, best_err = None, np.inf
    for perm in itertools.permutations([(0, 1), (1, 2), (2, 0)]):
        err = 0.0
        for k, (a, b) in enumerate(perm):
            mid = 0.5*(v[:, a] + v[:, b])
            err += np.abs(mids[:, k] - mid).max()
        if err < best_err:
            best, best_err = perm, err
    return best, best_err


def p2_eval(p, ue, x0, ys, order):
    """Evaluate u and du/dy at (x0, ys) inside one triangle."""
    P0, P1, P2 = p[0], p[1], p[2]
    Mm = np.array([[P1[0]-P0[0], P2[0]-P0[0]],
                   [P1[1]-P0[1], P2[1]-P0[1]]])
    try:
        Minv = np.linalg.inv(Mm)
    except np.linalg.LinAlgError:
        return None
    rel = np.column_stack([np.full_like(ys, x0) - P0[0], ys - P0[1]])
    L = rel @ Minv.T                      # (n,2): L1, L2
    L1 = L[:, 0]; L2 = L[:, 1]; L3 = 1.0 - L1 - L2
    inside = (L1 >= -1e-12) & (L2 >= -1e-12) & (L3 >= -1e-12)
    if not inside.any():
        return None
    dL1 = Minv[0, 1]; dL2 = Minv[1, 1]; dL3 = -(dL1 + dL2)
    Nc = {(0, 1): 4*L1*L2, (1, 2): 4*L2*L3, (2, 0): 4*L3*L1}
    dNc = {(0, 1): 4*(L1*dL2 + L2*dL1), (1, 2): 4*(L2*dL3 + L3*dL2),
           (2, 0): 4*(L3*dL1 + L1*dL3)}
    N = np.column_stack([L1*(2*L1-1), L2*(2*L2-1), L3*(2*L3-1),
                         Nc[order[0]], Nc[order[1]], Nc[order[2]]])
    dN = np.column_stack([(4*L1-1)*dL1, (4*L2-1)*dL2, (4*L3-1)*dL3,
                          dNc[order[0]], dNc[order[1]], dNc[order[2]]])
    return (N @ ue)[inside], (dN @ ue)[inside], ys[inside]


def profile(nodes, tris, u, x0, order, npts=41):
    p_all = nodes[tris]
    xmin = p_all[:, :, 0].min(1); xmax = p_all[:, :, 0].max(1)
    sel = np.flatnonzero((xmin <= x0) & (xmax >= x0))
    ys_all, us_all, ds_all = [], [], []
    for e in sel:
        p = p_all[e, :3]
        ylo = p[:, 1].min(); yhi = p[:, 1].max()
        ys = np.linspace(ylo, yhi, npts)
        out = p2_eval(p, u[tris[e]], x0, ys, order)
        if out is None:
            continue
        uu, dd, yy = out
        ys_all.append(yy); us_all.append(uu); ds_all.append(dd)
    if not ys_all:
        return None
    ys = np.concatenate(ys_all); us = np.concatenate(us_all); ds = np.concatenate(ds_all)
    o = np.argsort(ys)
    ys, us, ds = ys[o], us[o], ds[o]
    keep = np.concatenate([[True], np.diff(ys) > 1e-12])
    return ys[keep], us[keep], ds[keep]


def roots(ys, us, eps):
    out = []
    f = us - eps
    s = np.sign(f)
    idx = np.flatnonzero(s[:-1]*s[1:] < 0)
    for i in idx:
        t = f[i]/(f[i] - f[i+1])
        out.append(float(ys[i] + t*(ys[i+1]-ys[i])))
    return out


def intervals_from_roots(r, tol=5e-3):
    if not r:
        return []
    r = sorted(r); seg = [[r[0]]]
    for v in r[1:]:
        if v - seg[-1][-1] > tol:
            seg.append([v])
        else:
            seg[-1].append(v)
    return [(s[0], s[-1]) for s in seg if len(s) >= 2]


def grad_at(ys, us, ds, y0, eps, win=0.004):
    m = np.abs(ys - y0) <= win
    if not m.any():
        return 0.0, 0.0
    # slope at the crossing: use the P2 derivative nearest the root
    j = np.argmin(np.abs(ys - y0))
    return float(ds[j]), float(us[j])


def main():
    order_info = {}
    data = {}
    for name, f in MESHES:
        z = np.load(f)
        nodes, tris, u = z["nodes"], z["tris"], z["u"]
        order, err = detect_order(nodes, tris)
        order_info[name] = dict(order=order, max_mid_err=float(err))
        data[name] = dict(nodes=nodes, tris=tris, u=u, order=order)

    print("mid-side ordering (max |node - edge midpoint|):")
    for k, v in order_info.items():
        print(f"  {k:4s} {v['order']}  err={v['max_mid_err']:.2e}")

    out = {"order": order_info, "cuts": {}}
    for x0 in CUTS:
        print(f"\n================ x = {x0} ================")
        rec = {}
        prof = {}
        for name, _ in MESHES:
            d = data[name]
            pr = profile(d["nodes"], d["tris"], d["u"], x0, d["order"])
            prof[name] = pr
            z = np.load(dict(MESHES)[name])
            iv = cut_intervals(z["nodes"], z["tris"], z["mat"], z["u"], x0, 1e-8, True)
            main, short = clean(iv)
            per_eps = {}
            for eps in EPSS:
                iv2 = cut_intervals(z["nodes"], z["tris"], z["mat"], z["u"], x0, eps, True)
                m2, s2 = clean(iv2)
                per_eps[f"{eps:g}"] = [[float(a), float(b)] for a, b in m2]
            rec[name] = dict(extract_1e8=[[float(a), float(b)] for a, b in main],
                             poly=per_eps, short=int(len(short)))
            print(f"\n{name}: extractor eps=1e-8 {[(round(a,4),round(b,4)) for a,b in main]}")
            for eps in EPSS[1:]:
                print(f"   eps={eps:g}: {[(round(a,4),round(b,4)) for a,b in per_eps[f'{eps:g}']]}")

        # endpoints, gradients, conditioning
        print(f"\n-- endpoint conditioning (eps = 1e-8) --")
        print(f"{'mesh':5s} {'endpoint':>9s} {'u(y*)':>11s} {'du/dy':>11s} {'dy per 1e-8':>13s}")
        for name, _ in MESHES:
            pr = prof[name]
            ys, us, ds = pr
            for (a, b) in rec[name]["poly"]["1e-08"]:
                for y0 in (a, b):
                    uu = float(np.interp(y0, ys, us))
                    g = float(np.interp(y0, ys, ds))
                    amb = 1e-8/abs(g) if abs(g) > 1e-12 else float("inf")
                    print(f"{name:5s} {y0:+9.4f} {uu:11.3e} {g:+11.3e} {amb:13.3e}")
                    rec[name].setdefault("cond", []).append(
                        dict(y=y0, u=uu, dudy=g, dy_for_1e8=amb))

        # gap maxima
        print(f"\n-- gap maxima --")
        for name, _ in MESHES:
            ivs = rec[name]["poly"]["1e-08"]
            if len(ivs) < 2:
                print(f"{name:5s} single interval")
                continue
            pr = prof[name]
            gmax = []
            for (a1, b1), (a2, b2) in zip(ivs[:-1], ivs[1:]):
                m = (pr[0] > b1) & (pr[0] < a2)
                gmax.append(float(pr[1][m].max()) if m.any() else float("nan"))
            rec[name]["gap_max"] = gmax
            print(f"{name:5s} " + "  ".join(f"{v:.2e}" for v in gmax))

        # cross-mesh differences on a common grid
        print(f"\n-- cross-mesh |du| on common y grid --")
        interp = {}
        for name, _ in MESHES:
            pr = prof[name]
            y, uu = pr[0], pr[1]
            m = (y >= YG[0]) & (y <= YG[-1])
            interp[name] = np.interp(YG, y[m], uu[m], left=np.nan, right=np.nan)
        pairs = [("L1b", "L1"), ("mid", "L1b"), ("mid", "L1"), ("r4", "mid")]
        cross = {}
        for a, b in pairs:
            d = np.abs(interp[a] - interp[b])
            if np.all(np.isnan(d)):
                continue
            j = int(np.nanargmax(d))
            cross[f"{a}-{b}"] = dict(max=float(np.nanmax(d)), y=float(YG[j]),
                                     rms=float(np.sqrt(np.nanmean(d**2))))
            print(f"  {a:4s}-{b:4s} max |du| = {cross[f'{a}-{b}']['max']:.3e} at y={YG[j]:+.4f}  "
                  f"rms = {cross[f'{a}-{b}']['rms']:.3e}")
        rec["cross"] = cross
        out["cuts"][f"{x0}"] = rec

    json.dump(out, open("p2_front_conditioning.json", "w"), indent=2)
    print("\nwrote p2_front_conditioning.json")


if __name__ == "__main__":
    main()
