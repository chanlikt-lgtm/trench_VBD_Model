#!/usr/bin/env python3
"""Full-device potential maps in the reference style, for L1 and L1b."""
from __future__ import annotations
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D

import p2_audit as A
import mesh_stage2b as M2

M = A.M
V_D = A.V_DRAIN

xs = np.linspace(0.0, 2.42, 400)
a = M2.a_of(xs)


def bounds(ax, color="w", lw=1.4, ls="--"):
    for s in (1, -1):
        ax.plot(xs, s*a, color=color, lw=lw, ls=ls)
        ax.plot(xs, s*(a + 0.10), color=color, lw=lw*0.7, ls=ls, alpha=0.8)
    ax.set_xlim(0, 2.42); ax.set_ylim(-0.675, 0.675)
    ax.set_aspect("equal"); ax.set_xlabel("x (um)"); ax.set_ylabel("y (um)")


def labels(ax):
    kw = dict(ha="center", va="center", color="w", fontsize=10,
              path_effects=[pe.withStroke(linewidth=2.5, foreground="k")])
    ax.text(1.6, 0.0, "Si (n-drift)", **kw)
    ax.text(1.6, 0.605, "SiO2", fontsize=9, **{k: v for k, v in kw.items() if k != "fontsize"})
    ax.text(0.32, 0.55, "G", fontsize=9, **{k: v for k, v in kw.items() if k != "fontsize"})
    ax.text(1.35, 0.605, "SG", fontsize=9, **{k: v for k, v in kw.items() if k != "fontsize"})


def load(f):
    z = np.load(f)
    u = z["u"]
    phi = z["phi"] if "phi" in z else V_D - u
    return z["nodes"], z["tris"][:, :3], z["mat"], phi, u


def phi_map(lab, f, out):
    nodes, tv, mat, phi, u = load(f)
    x, y = nodes[:, 0], nodes[:, 1]
    fig, ax = plt.subplots(figsize=(10, 5.2))
    tcf = ax.tricontourf(x, y, tv, phi, levels=40, cmap="viridis")
    bounds(ax, color="w"); labels(ax)
    fig.colorbar(tcf, ax=ax, label="phi (V)")
    ax.set_title(f"{lab}: potential")
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)
    print(out)


def neutral_map(lab, f, out):
    nodes, tv, mat, phi, u = load(f)
    x, y = nodes[:, 0], nodes[:, 1]
    fig, ax = plt.subplots(figsize=(10, 5.2))
    simask = mat == M.SI
    face = np.where(simask, 0.82, 0.97)
    ax.tripcolor(x, y, tv, facecolors=face, shading="flat", cmap="gray", vmin=0, vmax=1)
    levels = [1e-8, 1e-4, 1e-2, 1e-1, 1.0]
    cs = ax.tricontour(x, y, tv, u, levels=levels, colors="k", linewidths=0.8)
    ax.clabel(cs, fmt=lambda v: f"{v:g}", fontsize=7)
    ax.tricontour(x, y, tv, u, levels=[1e-8], colors="r", linewidths=2.2)
    bounds(ax, color="k", lw=1.2); labels(ax)
    handles = [Line2D([0], [0], color="r", lw=2.2, label="u = 1e-8 (front)"),
               Line2D([0], [0], color="k", lw=0.8, label="u contours")]
    ax.legend(handles=handles, loc="upper right", fontsize=8)
    ax.set_title(f"{lab}: neutral region and depletion front")
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)
    print(out)


def main():
    phi_map("P2 L1b (h_loc 0.0055)", "p2_L1b_scaled.npz", "p2_L1b_phi.png")
    phi_map("P2 L1 (h_loc 0.0040)", "p2_L1_scaled.npz", "p2_L1_phi.png")
    neutral_map("P2 L1b", "p2_L1b_scaled.npz", "p2_L1b_neutral.png")
    neutral_map("P2 L1", "p2_L1_scaled.npz", "p2_L1_neutral.png")


if __name__ == "__main__":
    main()
