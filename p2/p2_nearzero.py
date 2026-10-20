#!/usr/bin/env python3
"""P2.1 near-zero-flux classification and pair-count reconciliation.

Part 0 reconciles the reported 810/170 against 896 total. That number came from
a single ad-hoc print in the previous session, where a *fraction* (0.8103) was
labelled "same-sign frac". This script recomputes both the fraction and the
integer counts from the same matched pairs and prints them so the two cannot be
confused again.

Part 1 classifies each matched quadrature point by
    D_scale = (|D_si| + |D_ox|) / 2
against a documented reference scale, at thresholds 1e-4, 1e-3, 1e-2.
Near-zero pairs are excluded only from the pointwise relative-mismatch
distribution; every sample stays in the integrated R_abs.
"""
from __future__ import annotations
import json
import numpy as np

import p2_interface as PI

THRESHOLDS = (1e-4, 1e-3, 1e-2)


def collect(flip=False):
    z = np.load("p2_r3_solution.npz")
    nodes, tris, mat, phi = z["nodes"], z["tris"], z["mat"], z["phi"]
    return PI.interface_flux(nodes, tris, mat, phi, flip_normal=flip)


def R_abs(recs):
    w = np.array([q["ds"] for q in recs])
    ds = np.array([q["d_si"] for q in recs])
    do = np.array([q["d_ox"] for q in recs])
    return float(np.sum(w * np.abs(ds - do)) /
                 np.sum(w * (np.abs(ds) + np.abs(do)) / 2))


def R_signed(recs):
    w = np.array([q["ds"] for q in recs])
    ds = np.array([q["d_si"] for q in recs])
    do = np.array([q["d_ox"] for q in recs])
    return float(abs(np.sum(w * (ds - do))) /
                 np.sum(w * (np.abs(ds) + np.abs(do)) / 2))


def main():
    recs = collect()
    ds = np.array([q["d_si"] for q in recs])
    do = np.array([q["d_ox"] for q in recs])
    w = np.array([q["ds"] for q in recs])
    dscale = (np.abs(ds) + np.abs(do)) / 2
    same = np.sign(ds) == np.sign(do)
    opp = ~same

    out = {}

    # ---------- reconciliation ----------
    wall = np.array([q["region"] == "wall" for q in recs])
    arc = ~wall
    out["counts"] = dict(
        total=int(len(recs)),
        wall=int(wall.sum()), arc=int(arc.sum()),
        same_sign_all=int(same.sum()), opp_sign_all=int(opp.sum()),
        same_plus_opp=int(same.sum() + opp.sum()),
        wall_same=int(same[wall].sum()), wall_opp=int(opp[wall].sum()),
        wall_same_plus_opp=int(same[wall].sum() + opp[wall].sum()),
        wall_same_fraction=float(same[wall].mean()),
        arc_same=int(same[arc].sum()), arc_opp=int(opp[arc].sum()),
        zero_pairs=int((ds == 0).sum() + (do == 0).sum()),
    )
    c = out["counts"]
    print("[counts] total={total} wall={wall} arc={arc}".format(**c))
    print(f"[counts] wall same={c['wall_same']} opp={c['wall_opp']} "
          f"sum={c['wall_same_plus_opp']} (fraction {c['wall_same_fraction']:.4f})")
    print(f"[counts] all  same={c['same_sign_all']} opp={c['opp_sign_all']} "
          f"sum={c['same_plus_opp']}")

    # ---------- reference scale ----------
    dref = float(np.median(dscale))
    out["reference"] = dict(d_ref=float(dref), rule="median of D_scale over all samples")
    print(f"[scale] D_ref (median D_scale) = {dref:.6f}")

    # ---------- opposite-sign minimum magnitude ----------
    out["opposite_sign"] = dict(
        count=int(opp.sum()),
        min_Dscale=float(dscale[opp].min()) if opp.any() else None,
        median_Dscale=float(np.median(dscale[opp])) if opp.any() else None,
        n_below_1e3=float((dscale[opp] < 1e-3*dref).sum()) if opp.any() else 0,
        wall=dict(count=int(opp[wall].sum()),
                  min_Dscale=float(dscale[wall & opp].min()) if (wall & opp).any() else None,
                  median_Dscale=float(np.median(dscale[wall & opp])) if (wall & opp).any() else None),
    )
    o = out["opposite_sign"]
    print(f"[opp] n={o['count']} min D_scale={o['min_Dscale']:.4e} "
          f"median={o['median_Dscale']:.4e}")

    # ---------- threshold sweep ----------
    rel = np.abs(ds - do) / np.maximum(dscale * 2, 1e-30)
    out["thresholds"] = {}
    for th in THRESHOLDS:
        cut = th * dref
        keep = dscale >= cut
        nz = ~keep
        entry = dict(
            threshold=th, cutoff=float(cut),
            near_zero=int(nz.sum()), kept=int(keep.sum()),
            near_zero_wall=int((nz & wall).sum()),
            near_zero_arc=int((nz & arc).sum()),
            rel_median=float(np.median(rel[keep])) if keep.any() else None,
            rel_mean=float(np.mean(rel[keep])) if keep.any() else None,
            rel_p95=float(np.percentile(rel[keep], 95)) if keep.any() else None,
            rel_max=float(np.max(rel[keep])) if keep.any() else None,
            R_abs_all=float(np.sum(w * np.abs(ds - do)) /
                            np.sum(w * (np.abs(ds) + np.abs(do)) / 2)),
        )
        if nz.any():
            entry["near_zero_x_range"] = [float(np.min([q["x"] for q, m in zip(recs, nz) if m])),
                                          float(np.max([q["x"] for q, m in zip(recs, nz) if m]))]
            entry["near_zero_y_range"] = [float(np.min([q["y"] for q, m in zip(recs, nz) if m])),
                                          float(np.max([q["y"] for q, m in zip(recs, nz) if m]))]
        out["thresholds"][str(th)] = entry
        print(f"[th {th:g}] cutoff={cut:.4e} near-zero={entry['near_zero']} "
              f"(wall {entry['near_zero_wall']}, arc {entry['near_zero_arc']}) "
              f"rel_med={entry['rel_median']:.4f} p95={entry['rel_p95']:.4f} "
              f"R_abs={entry['R_abs_all']:.6f}")

    out["integrated"] = dict(R_abs=R_abs(recs), R_signed=R_signed(recs),
                             wall_R_abs=R_abs([q for q in recs if q["region"] == "wall"]),
                             arc_R_abs=R_abs([q for q in recs if q["region"] == "arc"]),
                             wall_R_signed=R_signed([q for q in recs if q["region"] == "wall"]),
                             arc_R_signed=R_signed([q for q in recs if q["region"] == "arc"]))
    print(f"[integrated] R_abs={out['integrated']['R_abs']:.10f} "
          f"R_signed={out['integrated']['R_signed']:.6f}")

    with open("p2_nearzero.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print("wrote p2_nearzero.json")


if __name__ == "__main__":
    main()
