#!/usr/bin/env python3
"""Stage 2B.2 front gate — depletion-front geometry on the three existing meshes.

The neutral set is the active obstacle set on free silicon nodes. Drain nodes
(u = 0 by Dirichlet construction) are excluded before any measurement, so the
prescribed boundary cannot inflate the area.

Front reconstruction from the P1 solution:
  * neutral nodes are those where the complementarity constraint is active,
    equivalently u = V_D - phi <= tol_u with the node not fixed;
  * a depth cut at x_d takes every neutral node with |x - x_d| <= band, and the
    front is reported as the most negative and most positive y among them;
  * the same cut is repeated at three values of tol_u to expose sensitivity to
    the numerical zero threshold.
"""
from __future__ import annotations
import json
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

import stage2b2 as S

DEPTHS = [2.05, 2.10, 2.15, 2.20, 2.30, 2.35, 2.40]
BAND = 0.02
TOLS = [1e-9, 1e-6, 1e-3]
RES = ["r5", "fA", "fB"]


def components_of(nodes, tris, mask):
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return 0, np.zeros(0, int), idx
    local = {int(n): k for k, n in enumerate(idx)}
    rows, cols = [], []
    for e in range(len(tris)):
        a, b, c = tris[e]
        for u, v in ((a, b), (b, c), (c, a)):
            if u in local and v in local:
                rows += [local[int(u)], local[int(v)]]
                cols += [local[int(v)], local[int(u)]]
    A = sp.coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(idx.size, idx.size))
    ncomp, lab = connected_components(A.tocsr(), directed=False)
    return ncomp, lab, idx


def p1_neutral_area(nodes, tris, mats, u):
    """Exact area of {u <= 0} for a piecewise-linear u, per element.

    The zero level set of a P1 function is piecewise linear, so clipping each
    triangle against u = 0 gives the exact measure of the P1 neutral set. This
    is independent of how many nodes happen to fall inside the patch, which is
    the weakness of a nodal-area (lumped) estimator.
    """
    total = 0.0
    for e in range(len(tris)):
        if mats[e] != S.SI:
            continue
        p = nodes[tris[e]]
        val = u[tris[e]]
        ins = val <= 0.0
        if ins.all():
            total += 0.5 * abs((p[1, 0] - p[0, 0]) * (p[2, 1] - p[0, 1]) -
                               (p[1, 1] - p[0, 1]) * (p[2, 0] - p[0, 0]))
            continue
        if not ins.any():
            continue
        poly = []
        for k in range(3):
            a, b = k, (k + 1) % 3
            if ins[a]:
                poly.append(p[a])
            if ins[a] != ins[b]:
                t = val[a] / (val[a] - val[b])
                poly.append(p[a] + t * (p[b] - p[a]))
        if len(poly) < 3:
            continue
        q = np.array(poly)
        xq, yq = q[:, 0], q[:, 1]
        total += 0.5 * abs(np.dot(xq, np.roll(yq, -1)) - np.dot(yq, np.roll(xq, -1)))
    return float(total)


def analyse(tag):
    nodes, tris, mats = S.load(f"{tag}_mesh.npz")
    sol = S.solve_u(nodes, tris, mats)
    comp = sol["comp"]
    x, y = nodes[:, 0], nodes[:, 1]

    # drain Dirichlet nodes, excluded from every measurement
    si_n = sol["si_n"]
    drain_n = {i for i in si_n if x[i] >= S.BOTTOM_DEPTH - 1e-9}
    active = sol["active"].copy()
    active[list(drain_n)] = False

    n_act = int(active.sum())
    nodal_area = np.zeros(len(nodes))
    si = mats == S.SI
    np.add.at(nodal_area, tris[si].T,
              np.broadcast_to(sol["area"][si] / 3.0, (3, int(si.sum()))))
    area = float(nodal_area[active].sum())
    si_area = float(nodal_area[list(si_n)].sum())
    p1_area = p1_neutral_area(nodes, tris, mats, sol["u"])

    bbox = None
    if n_act:
        bbox = dict(xmin=float(x[active].min()), xmax=float(x[active].max()),
                    ymin=float(y[active].min()), ymax=float(y[active].max()))

    ncomp, lab, idx = components_of(nodes, tris, active)
    drain_conn = [False] * max(ncomp, 0)
    bridging_elements = 0
    bridging_active_nodes = set()
    for e in range(len(tris)):
        vs = tris[e]
        has_drain = any(int(v) in drain_n for v in vs)
        if not has_drain:
            continue
        acts = [int(v) for v in vs if active[v]]
        if acts:
            bridging_elements += 1
            for v in acts:
                bridging_active_nodes.add(v)
            for v in acts:
                k = int(np.where(idx == v)[0][0])
                drain_conn[int(lab[k])] = True
    gap_to_drain = float(S.BOTTOM_DEPTH - x[active].max()) if n_act else None

    fronts = {}
    for xd in DEPTHS:
        sel = std = None
        row = {"x": xd}
        for side, pick in (("neg", np.min), ("pos", np.max)):
            band = np.flatnonzero(active & (np.abs(x - xd) <= BAND))
            row[side] = float(pick(y[band])) if band.size else None
        fronts[f"{xd:.2f}"] = row

    sensitivity = {}
    for tol in TOLS:
        mask = active | ((sol["u"] <= tol) & np.array([i in si_n for i in range(len(nodes))])
                         & ~np.isin(np.arange(len(nodes)), list(drain_n)))
        sens = {}
        for xd in DEPTHS:
            band = np.flatnonzero(mask & (np.abs(x - xd) <= BAND))
            sens[f"{xd:.2f}"] = dict(neg=float(y[band].min()) if band.size else None,
                                     pos=float(y[band].max()) if band.size else None)
        sensitivity[f"tol={tol:g}"] = sens

    asym = {}
    for k, row in fronts.items():
        if row["neg"] is not None and row["pos"] is not None:
            asym[k] = float(abs(abs(row["neg"]) - abs(row["pos"])))

    return dict(
        tag=tag, triangles=int(len(tris)),
        active_nodes=n_act, active_nodes_on_drain_plane=int(
            sum(1 for i in np.flatnonzero(active) if x[i] >= S.BOTTOM_DEPTH - 1e-9)),
        neutral_area_um2=area, silicon_area_um2=si_area,
        p1_neutral_area_um2=p1_area,
        p1_fraction=p1_area / si_area if si_area else None,
        neutral_fraction=area / si_area if si_area else None,
        bbox=bbox,
        components=dict(count=int(ncomp), drain_connected=bool(any(drain_conn)),
                        bridging_elements=int(bridging_elements),
                        bridging_active_nodes=int(len(bridging_active_nodes)),
                        gap_to_drain_um=gap_to_drain,
                        sizes=sorted([int((lab == c).sum()) for c in range(ncomp)], reverse=True)),
        fronts=fronts,
        asymmetry_um=asym,
        threshold_sensitivity=sensitivity,
        solver=dict(iterations=comp["iterations"],
                    primal_min_u=comp["primal_min_u_si"],
                    dual_min_r=comp["dual_min_r_si"],
                    complementarity_norm=comp["complementarity_norm"]),
    )


def main():
    out = {t: analyse(t) for t in RES}
    with open("front_gate.json", "w") as f:
        json.dump(out, f, indent=2)

    for t in RES:
        r = out[t]
        b = r["bbox"]
        print(f"\n=== {t}: {r['triangles']} tri, {r['active_nodes']} active "
              f"({r['neutral_area_um2']:.5f} um2, {100*r['neutral_fraction']:.3f}% of Si) ===")
        print(f"  bbox  x[{b['xmin']:.4f},{b['xmax']:.4f}]  y[{b['ymin']:.4f},{b['ymax']:.4f}]"
              if b else "  bbox  N/A")
        print(f"  components {r['components']['count']} drain_connected="
              f"{r['components']['drain_connected']} sizes={r['components']['sizes']}")
        print(f"  active nodes on drain plane: {r['active_nodes_on_drain_plane']}")
        print("  front (neg / pos y):")
        for k, v in r["fronts"].items():
            n = "N/A" if v["neg"] is None else f"{v['neg']:+.4f}"
            p = "N/A" if v["pos"] is None else f"{v['pos']:+.4f}"
            print(f"    x={k}  {n}  {p}")
        print("  asymmetry:", {k: round(v, 4) for k, v in r["asymmetry_um"].items()})

    print("\nresolution dependence (neutral area, bbox, component count):")
    print(f"{'tag':>4} {'tri':>7} {'lump_um2':>10} {'lump%':>7} {'P1_um2':>9} {'P1%':>7} {'ncomp':>6} "
          f"{'xmin':>7} {'xmax':>7} {'ymin':>7} {'ymax':>7}")
    for t in RES:
        r = out[t]; b = r["bbox"] or dict(xmin=np.nan, xmax=np.nan, ymin=np.nan, ymax=np.nan)
        print(f"{t:>4} {r['triangles']:>7} {r['neutral_area_um2']:>10.5f} "
              f"{100*r['neutral_fraction']:>7.3f} {r['p1_neutral_area_um2']:>9.5f} "
              f"{100*r['p1_fraction']:>7.3f} {r['components']['count']:>6} "
              f"{b['xmin']:>7.4f} {b['xmax']:>7.4f} {b['ymin']:>7.4f} {b['ymax']:>7.4f}")


if __name__ == "__main__":
    main()
