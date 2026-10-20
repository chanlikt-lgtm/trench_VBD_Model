#!/usr/bin/env python3
"""P2.3 neutral-geometry audit on the polished solution.

- documented zero-set classes (fully neutral, any-single-zero, drain-adjacent)
- connected components of the fully neutral set via shared nodes
- drain connectivity
- front intersections at x = 2.35, 2.38, 2.39, 2.40, 2.41 (upper and lower y)
- transverse symmetry
- interface flux by region after activation
- comparison with the P1 archived active area
"""
from __future__ import annotations
import json
import sys
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

import p2_audit as A
import p2_interface as PI

M = A.M
V_D = A.V_DRAIN
TOL = 1e-8


def main():
    FN = sys.argv[1] if len(sys.argv) > 1 else "p2_p3_r3_polished.npz"
    OUT = sys.argv[2] if len(sys.argv) > 2 else "p2_p3_geometry.json"
    z = np.load(FN)
    nodes, tris, mat = z["nodes"], z["tris"], z["mat"]
    phi, u = z["phi"], z["u"]
    area = np.abs(A.element_geom(nodes, tris)[2])

    cxs = nodes[tris][:, :3, 0].mean(axis=1)
    si = mat == M.SI
    body = si & (cxs <= A.BODY_BOTTOM + 1e-12)
    charged = si & ~body
    ue = u[tris]

    fully = charged & (np.abs(ue) < TOL).all(axis=1)
    anyz = charged & (np.abs(ue) < TOL).any(axis=1)
    drain_node = nodes[:, 0] >= M.BOTTOM_DEPTH - 1e-9
    elem_drain = charged & drain_node[tris].any(axis=1)

    print(f"[classes] charged elements={int(charged.sum())} area={area[charged].sum():.9f}")
    print(f"  fully neutral : {int(fully.sum()):5d} elements  area={area[fully].sum():.9f}")
    print(f"  any zero node : {int(anyz.sum()):5d} elements  area={area[anyz].sum():.9f}")
    print(f"  drain-adjacent: {int(elem_drain.sum()):5d} elements  area={area[elem_drain].sum():.9f}")
    print(f"  fully neutral & drain-adjacent: "
          f"{int((fully & elem_drain).sum())} elements  "
          f"area={area[fully & elem_drain].sum():.9f}")
    print(f"  fully neutral & NOT drain-adjacent: "
          f"{int((fully & ~elem_drain).sum())} elements  "
          f"area={area[fully & ~elem_drain].sum():.9f}")

    # components among fully neutral elements via shared nodes
    els = np.flatnonzero(fully)
    local = {int(e): k for k, e in enumerate(els)}
    node2 = {}
    for e in els:
        for n in tris[e]:
            node2.setdefault(int(n), []).append(int(e))
    rows, cols = [], []
    for n, es in node2.items():
        for i in es:
            for j in es:
                if i != j:
                    rows.append(local[i]); cols.append(local[j])
    Adj = sp.coo_matrix((np.ones(len(rows)), (rows, cols)),
                        shape=(len(els), len(els))).tocsr()
    ncomp, lab = connected_components(Adj, directed=False) if len(els) else (0, np.array([]))
    print(f"\n[components] fully neutral components={ncomp}")
    for c in range(ncomp):
        m = lab == c
        es = els[m]
        a = area[es].sum()
        pts = nodes[np.unique(tris[es])]
        touches = bool(drain_node[np.unique(tris[es])].any())
        print(f"  comp {c}: elements={int(m.sum()):5d} area={a:.9f} "
              f"x[{pts[:,0].min():.4f},{pts[:,0].max():.4f}] "
              f"y[{pts[:,1].min():.4f},{pts[:,1].max():.4f}] drain={touches}")

    # drain connectivity fraction
    drain_si = np.unique(tris[si][drain_node[tris[si]]])
    if len(drain_si):
        frac = float(np.mean(np.abs(u[drain_si]) < TOL))
        print(f"[drain] Si drain-plane nodes={len(drain_si)}  |u|<tol fraction={frac:.4f}")

    # front cuts
    def eval_cut(x0, ys):
        cidx = np.flatnonzero(charged)
        Aa = nodes[tris[cidx, 0]]; Bb = nodes[tris[cidx, 1]]; Cc = nodes[tris[cidx, 2]]
        det = ((Bb[:, 0]-Aa[:, 0])*(Cc[:, 1]-Aa[:, 1]) -
               (Cc[:, 0]-Aa[:, 0])*(Bb[:, 1]-Aa[:, 1]))
        vals = []
        for y in ys:
            L1 = ((Bb[:, 1]-Cc[:, 1])*(x0-Cc[:, 0]) +
                  (Cc[:, 0]-Bb[:, 0])*(y-Cc[:, 1]))/det
            L2 = ((Cc[:, 1]-Aa[:, 1])*(x0-Cc[:, 0]) +
                  (Aa[:, 0]-Cc[:, 0])*(y-Cc[:, 1]))/det
            L3 = 1-L1-L2
            m = (L1 >= -1e-9) & (L2 >= -1e-9) & (L3 >= -1e-9)
            if not m.any():
                vals.append(np.nan); continue
            j = int(np.flatnonzero(m)[0])
            e = int(cidx[j])
            vals.append(float(np.dot(u[tris[e]], A.shape_vals(L1[j], L2[j]))))
        return np.array(vals)

    print("\n[front cuts]")
    ys = np.arange(-0.6, 0.6001, 0.002)
    out_cuts = {}
    for x0 in (2.35, 2.38, 2.39, 2.40, 2.41):
        v = eval_cut(x0, ys)
        ok = ~np.isnan(v)
        neu = ok & (np.abs(v) < TOL)
        # contiguous intervals
        intervals = []
        if neu.any():
            idx = np.flatnonzero(neu)
            start = idx[0]; prev = idx[0]
            for k in idx[1:]:
                if k != prev+1:
                    intervals.append((ys[start], ys[prev])); start = k
                prev = k
            intervals.append((ys[start], ys[prev]))
        out_cuts[x0] = intervals
        print(f"  x={x0:.2f}: Si extent y[{ys[ok].min():+.3f},{ys[ok].max():+.3f}] "
              f"neutral intervals={[(round(a,3), round(b,3)) for a, b in intervals]}")

    # transverse symmetry
    dmax = 0.0; dloc = None
    for x0 in np.arange(2.33, 2.421, 0.01):
        yy = np.arange(0.0, 0.55, 0.01)
        vp = eval_cut(x0, yy); vm = eval_cut(x0, -yy)
        m = ~np.isnan(vp) & ~np.isnan(vm)
        if m.any():
            d = np.abs(vp[m]-vm[m])
            k = int(np.argmax(d))
            if d[k] > dmax:
                dmax = float(d[k]); dloc = (float(x0), float(yy[m][k]))
    print(f"\n[symmetry] max |u(x,y)-u(x,-y)| = {dmax:.3e} V at {dloc}")

    # interface flux by region after activation
    recs = PI.interface_flux(nodes, tris, mat, phi)

    def st(rr):
        if not rr:
            return None
        w = np.array([q["ds"] for q in rr]); ds = np.array([q["d_si"] for q in rr])
        do = np.array([q["d_ox"] for q in rr])
        den = np.sum(w*(np.abs(ds)+np.abs(do))/2)
        return dict(n=len(rr), R_abs=float(np.sum(w*np.abs(ds-do))/den),
                    R_sgn=float(abs(np.sum(w*(ds-do)))/den))

    regions = {"wall": [q for q in recs if q["x"] <= 2.04],
               "arc": [q for q in recs if 2.04 < q["x"] <= 2.14],
               "cap": [q for q in recs if q["x"] > 2.14]}
    print("\n[interface after activation]")
    res = {}
    for k, rr in regions.items():
        s = st(rr); res[k] = s
        if s:
            print(f"  {k:4s} n={s['n']:4d} R_abs={s['R_abs']:.6f} R_sgn={s['R_sgn']:.6f}")

    out = dict(fully_area=float(area[fully].sum()), fully_elements=int(fully.sum()),
               any_area=float(area[anyz].sum()), any_elements=int(anyz.sum()),
               drain_adjacent_elements=int(elem_drain.sum()),
               components=int(ncomp), cuts={str(k): v for k, v in out_cuts.items()},
               symmetry_max=dmax, symmetry_loc=dloc, interface=res)
    json.dump(out, open(OUT, "w"), indent=2)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
