#!/usr/bin/env python3
"""Stage 2B.3 — front reconstruction from the P1 obstacle solution.

Two complementary reconstructions, as required:

  (1) EXACT discrete zero set.  The obstacle enforces u_i >= 0, so a conventional
      marching-triangles sign test has nothing to bracket.  The exact zero set is
      the union of triangles whose three nodes all have u = 0 (up to tol), plus
      their connecting edges.  Its boundary is traced as a polyline.

  (2) SMALL POSITIVE LEVEL SETS  u = delta, delta -> 0+.  These are offsets from
      the true front, not the front, and are reported as such.

Front positions are line intersections of the boundary with x = const, not the
largest neutral-node coordinate in a window.
"""
from __future__ import annotations
import json
import numpy as np

import stage2b2 as S

CASES = ["r5", "fA", "fB"]
CUTS = [2.35, 2.38, 2.39, 2.40, 2.41]
DELTAS = [1e-12, 1e-6, 1e-4, 1e-3, 1e-2]
ZERO_TOL = 1e-9


def boundary_edges(nodes, tris, mats, u):
    """Ordered boundary polyline(s) of the exact zero set."""
    zero = (mats == S.SI) & (u[tris] <= ZERO_TOL).all(axis=1)
    zt = tris[zero]
    if len(zt) == 0:
        return [], 0

    edge_count = {}
    for e in zt:
        for a, b in ((e[0], e[1]), (e[1], e[2]), (e[2], e[0])):
            k = (int(min(a, b)), int(max(a, b)))
            edge_count[k] = edge_count.get(k, 0) + 1

    bedges = [k for k, c in edge_count.items() if c == 1]
    adj = {}
    for a, b in bedges:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)

    loops = []
    seen = set()
    for start in adj:
        if start in seen:
            continue
        loop = [start]
        seen.add(start)
        prev, cur = None, start
        while True:
            nxt = [n for n in adj[cur] if n != prev]
            if not nxt:
                break
            prev, cur = cur, nxt[0]
            if cur == start:
                break
            if cur in seen:
                break
            loop.append(cur)
            seen.add(cur)
        if len(loop) >= 3:
            loops.append(loop)
    return loops, int(zero.sum())


def cut_of_loop(nodes, loop, xd):
    """All y where the closed loop crosses x = xd, split into negative/positive."""
    ys = []
    n = len(loop)
    for k in range(n):
        p, q = nodes[loop[k]], nodes[loop[(k + 1) % n]]
        if (p[0] - xd) * (q[0] - xd) < 0:
            t = (xd - p[0]) / (q[0] - p[0])
            ys.append(p[1] + t * (q[1] - p[1]))
    neg = [v for v in ys if v < 0]
    pos = [v for v in ys if v >= 0]
    return (max(neg) if neg else None), (min(pos) if pos else None), len(ys)


def level_cut(nodes, tris, u, xd, delta):
    """Intersections of u = delta with the line x = xd (sub-element positions)."""
    out = []
    for e in range(len(tris)):
        p = nodes[tris[e]]
        val = u[tris[e]] - delta
        ins = val >= 0
        if ins.all() or not ins.any():
            continue
        pts = []
        for k in range(3):
            a, b = k, (k + 1) % 3
            if ins[a] == ins[b]:
                continue
            t = val[a] / (val[a] - val[b])
            pts.append(p[a] + t * (p[b] - p[a]))
        for i in range(len(pts) - 1):
            P, Q = pts[i], pts[i + 1]
            if (P[0] - xd) * (Q[0] - xd) <= 0 and abs(Q[0] - P[0]) > 1e-15:
                t = (xd - P[0]) / (Q[0] - P[0])
                out.append(P[1] + t * (Q[1] - P[1]))
    neg = [v for v in out if v < 0]
    pos = [v for v in out if v >= 0]
    return (max(neg) if neg else None), (min(pos) if pos else None), len(out)


def main():
    res = {}
    for tag in CASES:
        nodes, tris, mats = S.load(f"{tag}_mesh.npz")
        sol = S.solve_u(nodes, tris, mats)
        u = sol["u"]
        si_area = sol["comp"]["si_area_um2"]

        loops, nzt = boundary_edges(nodes, tris, mats, u)
        main_loop = max(loops, key=len) if loops else []
        ncomp = len(loops)

        exact = {}
        for xd in CUTS:
            neg, pos, nint = cut_of_loop(nodes, main_loop, xd)
            exact[f"{xd:.2f}"] = dict(neg=neg, pos=pos, intersections=nint)

        sens = {}
        for d in DELTAS:
            sens[f"{d:g}"] = {f"{xd:.2f}": level_cut(nodes, tris, u, xd, d)
                              for xd in CUTS}

        ys = nodes[main_loop, 1] if main_loop else np.array([])
        res[tag] = dict(
            triangles=int(len(tris)),
            zero_triangles=nzt,
            boundary_loops=ncomp,
            max_transverse_extent=float(np.max(np.abs(ys))) if ys.size else None,
            cut_positions=exact,
            level_sensitivity={k: {x: v for x, v in d.items()} for k, d in sens.items()},
            solver=dict(primal=sol["comp"]["primal_min_u_si"],
                        dual=sol["comp"]["dual_min_r_si"],
                        comp=sol["comp"]["complementarity_norm"]),
        )

    with open("front_reconstruct.json", "w") as f:
        json.dump(res, f, indent=2)

    print(f"{'cut':>6} " + " ".join(f"{t:>22}" for t in CASES))
    print("-" * 78)
    for xd in CUTS:
        row = f"x={xd:>4.2f} "
        for t in CASES:
            d = res[t]["cut_positions"][f"{xd:.2f}"]
            neg = "N/A" if d["neg"] is None else f"{d['neg']:+.4f}"
            pos = "N/A" if d["pos"] is None else f"{d['pos']:+.4f}"
            row += f" {neg}/{pos} ({d['intersections']})".ljust(23)
        print(row)

    print("\nzero triangles / loops / max |y| extent:")
    for t in CASES:
        r = res[t]
        print(f"  {t}: {r['zero_triangles']} tri  loops={r['boundary_loops']}  "
              f"max|y|={r['max_transverse_extent']:.5f}")

    print("\nlevel-set sensitivity at x=2.40 (neg/pos as delta -> 0):")
    for t in CASES:
        print(f"  {t}:", end="")
        for d in DELTAS:
            neg, pos, _ = res[t]["level_sensitivity"][f"{d:g}"]["2.40"]
            print(f"  d={d:g}:{neg:+.4f}/{pos:+.4f}" if neg is not None
                  else f"  d={d:g}:N/A", end="")
        print()


if __name__ == "__main__":
    main()
