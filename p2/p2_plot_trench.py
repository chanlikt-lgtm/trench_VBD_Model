#!/usr/bin/env python3
"""Plot the solved potential in the trench window for L1 and L1b."""
from __future__ import annotations
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

import p2_audit as A

V_D = A.V_DRAIN
X0, X1, Y0, Y1 = 2.08, 2.56, -0.34, 0.34
EPS = 1e-8


def load(f):
    z = np.load(f)
    nodes, tris, u = z["nodes"], z["tris"], z["u"]
    phi = V_D - u
    cx = nodes[tris][:, :3, 0].mean(1)
    cy = nodes[tris][:, :3, 1].mean(1)
    keep = (cx >= X0) & (cx <= X1) & (cy >= Y0) & (cy <= Y1)
    t = tris[keep]
    used = np.unique(t.ravel())
    o2n = np.full(len(nodes), -1, dtype=int)
    o2n[used] = np.arange(len(used))
    return nodes[used], o2n[t][:, :3], phi[used], u[used]


def fig_maps():
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 7.6), sharex=True, sharey=True,
                             constrained_layout=True)
    for ax, (lab, f) in zip(axes, (("L1b, h_loc 0.0055", "p2_L1b_scaled.npz"),
                                   ("L1, h_loc 0.0040", "p2_L1_scaled.npz"))):
        n, t, phi, u = load(f)
        tri = mtri.Triangulation(n[:, 0], n[:, 1], t)
        cs = ax.tricontourf(tri, phi, levels=40, cmap="viridis")
        ax.tricontour(tri, u, levels=[EPS], colors="w", linewidths=0.9)
        ax.set_aspect("equal", adjustable="box", anchor="C")
        ax.set_title(f"φ in the trench window — {lab}", fontsize=9.5)
        ax.set_ylabel("y (µm)")
        cb = fig.colorbar(cs, ax=ax, fraction=0.035, pad=0.02)
        cb.set_label("φ (V)", fontsize=8)
        cb.ax.tick_params(labelsize=7)
    axes[-1].set_xlabel("x (µm)")
    fig.savefig("fig6_trench_phi.png")
    print("fig6_trench_phi.png")


def fig_boundary():
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    for lab, f, c in (("L1b", "p2_L1b_scaled.npz", "#1f6fb4"),
                      ("L1", "p2_L1_scaled.npz", "#c0392b")):
        n, t, phi, u = load(f)
        tri = mtri.Triangulation(n[:, 0], n[:, 1], t)
        ax.tricontour(tri, u, levels=[EPS], colors=c, linewidths=1.3)
        ax.plot([], [], color=c, lw=1.3, label=f"{lab} neutral boundary (u = 1e-8 V)")
    for xc in (2.35, 2.41):
        ax.axvline(xc, color="0.4", ls="--", lw=0.8)
        ax.annotate(f"x = {xc}", (xc, Y1), textcoords="offset points",
                    xytext=(3, -10), fontsize=8, color="0.3")
    ax.set_xlim(X0, X1); ax.set_ylim(Y0, Y1)
    ax.set_aspect("equal")
    ax.set_xlabel("x (µm)"); ax.set_ylabel("y (µm)")
    ax.set_title("Neutral-region boundary in the trench — L1b vs L1", fontsize=10)
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig("fig7_trench_boundary.png")
    print("fig7_trench_boundary.png")


def fig_u_profile():
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.9), sharey=True)
    for ax, xc in zip(axes, (2.35, 2.41)):
        for lab, f, c in (("L1b", "p2_L1b_scaled.npz", "#1f6fb4"),
                          ("L1", "p2_L1_scaled.npz", "#c0392b")):
            z = np.load(f)
            n, t, u = z["nodes"], z["tris"], z["u"]
            pts = []
            for e in t:
                p = n[e][:, :3]
                if p[:, 0].min() <= xc <= p[:, 0].max():
                    ymid = p[:, 1].mean()
                    if Y0 <= ymid <= Y1:
                        uu = u[e][:3].mean()
                        pts.append((ymid, uu))
            if pts:
                arr = np.array(sorted(pts))
                ax.plot(arr[:, 1], arr[:, 0], ".", ms=2.2, color=c, alpha=0.7, label=lab)
        ax.axvline(EPS, color="k", ls=":", lw=0.8)
        ax.set_xscale("symlog", linthresh=1e-8)
        ax.set_xlabel("u = V_D − φ  (V, symlog)")
        ax.set_title(f"x = {xc} µm", fontsize=9.5)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("y (µm)")
    axes[0].set_ylim(Y0, Y1)
    axes[0].legend(fontsize=7.5, loc="lower left", framealpha=0.85)
    fig.suptitle("Obstacle variable across the cut (element samples)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig("fig8_u_profiles.png")
    print("fig8_u_profiles.png")


if __name__ == "__main__":
    fig_maps()
    fig_boundary()
    fig_u_profile()
