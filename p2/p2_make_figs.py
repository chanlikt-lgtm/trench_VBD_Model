#!/usr/bin/env python3
"""Generate the P2.3 comparison figures from the certified solutions."""
from __future__ import annotations
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from p2_exact_metrics import cut_intervals
from p2_topology_clean import clean

plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3,
                     "figure.dpi": 160})

# mesh label, solution file, realised h in window A core (median), area interval
MESHES = [
    ("r3",  "p2_p3_r3_polished.npz",    0.01264, (0.020926837673, 0.020930251232)),
    ("r4",  "p2_p3_r4_polished.npz",    0.00909, (0.021574089791, 0.021578224392)),
    ("mid", "p2_p3_rmid_polished.npz",  0.00780, (0.021789286073, 0.021793581617)),
    ("L1b", "p2_L1b_scaled.npz",        0.00524, (0.021842101270, 0.021845749070)),
    ("L1",  "p2_L1_scaled.npz",         0.00399, (0.022071011922, 0.022075109085)),
]
COL = {"r3": "0.55", "r4": "0.35", "mid": "0.15", "L1b": "#1f6fb4", "L1": "#c0392b"}


def cuts(x0, eps):
    out = {}
    for name, f, h, _ in MESHES:
        z = np.load(f)
        iv = cut_intervals(z["nodes"], z["tris"], z["mat"], z["u"], x0, eps, True)
        main, short = clean(iv)
        out[name] = sorted([tuple(map(float, t)) for t in main])
    return out


def fig_area():
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    for name, _, h, (lo, up) in MESHES:
        ax.errorbar(h, 0.5*(lo+up), yerr=[[0.5*(up-lo)], [0.5*(up-lo)]],
                    fmt="o", ms=7, color=COL[name], capsize=4, lw=1.2)
        ax.annotate(name, (h, 0.5*(lo+up)), textcoords="offset points",
                    xytext=(8, -10), fontsize=9, color=COL[name])
    ax.set_xscale("log")
    ax.invert_xaxis()
    ax.margins(x=0.18)
    ax.set_xlabel("realised element size in window A (µm, median)  —  finer →")
    ax.set_ylabel("bounded neutral area (µm²)")
    ax.set_title("P2.3 neutral area vs realised mesh size\n"
                 "error bars are the Bernstein enclosure width", fontsize=10)
    fig.tight_layout()
    fig.savefig("fig1_area_vs_h.png")
    print("fig1_area_vs_h.png")


def fig_increments():
    names = [m[0] for m in MESHES]
    lab, lo, up = [], [], []
    for i in range(1, len(MESHES)):
        a = MESHES[i-1][3]; b = MESHES[i][3]
        lab.append(f"{names[i-1]}→{names[i]}")
        lo.append(b[0]-a[1]); up.append(b[1]-a[0])
    lo = np.array(lo); up = np.array(up)
    mid = 0.5*(lo+up); err = 0.5*(up-lo)
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.bar(range(len(lab)), mid, yerr=err, capsize=5,
           color=["0.6", "0.45", "0.3", "#1f6fb4"], edgecolor="k", lw=0.6)
    for i, (m, e) in enumerate(zip(mid, err)):
        ax.annotate(f"{m:.2e}\n(±{e:.1e})", (i, m), textcoords="offset points",
                    xytext=(0, 6), ha="center", fontsize=8)
    ax.set_xticks(range(len(lab)))
    ax.set_xticklabels(lab)
    ax.set_ylabel("neutral-area increment (µm²)")
    ax.set_title("Interval-bounded area increments\n"
                 "increment reduction reverses at the finest local refinement", fontsize=10)
    ax.set_ylim(0, max(up)*1.28)
    fig.tight_layout()
    fig.savefig("fig2_area_increments.png")
    print("fig2_area_increments.png")


def fig_front(x0, fname, title):
    data = cuts(x0, 1e-8)
    names = [m[0] for m in MESHES]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for i, name in enumerate(names):
        spans = [(a, b - a) for (a, b) in data[name]]
        ax.broken_barh(spans, (i - 0.28, 0.56), facecolors=COL[name],
                       edgecolors="k", linewidth=0.6, alpha=0.55)
        lo = min(a for a, b in data[name]); hi = max(b for a, b in data[name])
        ax.plot([i, i], [lo, hi], color=COL[name], lw=1.0, ls=":")
        ax.annotate(f"{len(spans)} interval{'s' if len(spans) > 1 else ''}",
                    (i, hi), textcoords="offset points", xytext=(0, 5),
                    ha="center", fontsize=7.5)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels([f"{n}\nh={m[2]:.4f}" for n, m in zip(names, MESHES)])
    ax.set_xlim(-0.6, len(names) - 0.4)
    ax.set_ylabel("y at the cut (µm)")
    ax.set_title(title, fontsize=10)
    ax.axhline(0, color="k", lw=0.6, ls=":")
    ax.plot([], [], "s", color="0.45", alpha=0.55, label="neutral interval (u < ε)")
    ax.legend(fontsize=7.5, loc="lower right")
    fig.tight_layout()
    fig.savefig(fname)
    print(fname)


def fig_edges():
    d35 = cuts(2.35, 1e-8); d41 = cuts(2.41, 1e-8)
    names = [m[0] for m in MESHES]
    hs = [m[2] for m in MESHES]
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.8), sharey=False)
    for ax, d, lab in ((axes[0], d35, "x = 2.35 µm"), (axes[1], d41, "x = 2.41 µm")):
        lower = [min(a for a, b in d[n]) for n in names]
        upper = [max(b for a, b in d[n]) for n in names]
        ax.plot(hs, lower, "o-", color="#1f6fb4", label="lower edge")
        ax.plot(hs, upper, "s--", color="#c0392b", label="upper edge")
        ax.set_xscale("log")
        ax.invert_xaxis()
        ax.set_xticks(hs)
        ax.set_xticklabels([f"{n}\n{h:.4f}" for n, h in zip(names, hs)], fontsize=7.5)
        ax.margins(x=0.12)
        ax.set_xlabel("realised h (µm) — finer →")
        ax.set_ylabel("y (µm)")
        ax.set_title(lab, fontsize=9.5)
        ax.legend(fontsize=7.5, loc="best")
    fig.suptitle("Neutral-front edges vs realised mesh size (ε = 1e-8 V)", fontsize=10)
    fig.tight_layout()
    fig.savefig("fig5_front_edges.png")
    print("fig5_front_edges.png")


if __name__ == "__main__":
    fig_area()
    fig_increments()
    fig_front(2.35, "fig3_front_x235.png", "Neutral intervals at x = 2.35 µm (ε = 1e-8 V)")
    fig_front(2.41, "fig4_front_x241.png", "Neutral intervals at x = 2.41 µm (ε = 1e-8 V)")
    fig_edges()
