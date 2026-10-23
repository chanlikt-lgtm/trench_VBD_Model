#!/usr/bin/env python3
"""|E| maps for a certified P2 solution, in the reference style.

E is evaluated at element centroids from the P2 gradient of phi = V_D - u.
Centroid values, not nodal maxima: the trench-tip maximum is not a converged
quantity and must not be quoted as one.
"""
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
SOL = "p2_B2X_scaled.npz"
XS = np.linspace(0.0, 2.42, 400)
AOF = M2.a_of(XS)


def load():
    z = np.load(SOL)
    nodes, tris, mat, u = z["nodes"], z["tris"], z["mat"], z["u"]
    phi = V_D - u
    _, det, area, g1, g2, g3 = A.element_geom(nodes, tris)
    L = 1.0/3.0
    G = np.zeros((len(tris), 6, 2))
    G[:, 0] = (4*L-1)*g1
    G[:, 1] = (4*L-1)*g2
    G[:, 2] = (4*L-1)*g3
    G[:, 3] = 4*(L*g2 + L*g1)
    G[:, 4] = 4*(L*g3 + L*g2)
    G[:, 5] = 4*(L*g1 + L*g3)
    E = -np.einsum('eip,ei->ep', G, phi[tris])
    Emv = np.linalg.norm(E, axis=1)/100.0        # V/um -> MV/cm
    return nodes, tris, mat, u, Emv


def bounds(ax, color="w", lw=1.4, ls="--"):
    for s in (1, -1):
        ax.plot(XS, s*AOF, color=color, lw=lw, ls=ls)
        ax.plot(XS, s*(AOF + 0.10), color=color, lw=lw*0.7, ls=ls, alpha=0.8)
    ax.set_xlim(0, 2.42); ax.set_ylim(-0.675, 0.675)
    ax.set_aspect("equal"); ax.set_xlabel("x (um)"); ax.set_ylabel("y (um)")


def labels(ax):
    kw = dict(ha="center", va="center", color="w", fontsize=10,
              path_effects=[pe.withStroke(linewidth=2.5, foreground="k")])
    ax.text(1.6, 0.0, "Si (n-drift)", **kw)
    ax.text(1.6, 0.605, "SiO2", fontsize=9, **{k: v for k, v in kw.items() if k != "fontsize"})
    ax.text(0.32, 0.55, "G", fontsize=9, **{k: v for k, v in kw.items() if k != "fontsize"})
    ax.text(1.35, 0.605, "SG", fontsize=9, **{k: v for k, v in kw.items() if k != "fontsize"})


def main():
    nodes, tris, mat, u, Emv = load()
    x, y = nodes[:, 0], nodes[:, 1]
    tv = tris[:, :3]
    print(f"centroid |E| max = {Emv.max():.3f} MV/cm "
          f"(centroid value; not a converged tip maximum)")

    for vmax, name in ((5.5, "figE1_full_0_5p5.png"), (2.0, "figE2_full_0_2p0.png")):
        fig, ax = plt.subplots(figsize=(10, 5.2))
        tp = ax.tripcolor(x, y, tv, facecolors=Emv, shading="flat",
                          cmap="turbo", vmin=0, vmax=vmax)
        bounds(ax)
        labels(ax)
        fig.colorbar(tp, ax=ax, label="|E| (MV/cm)")
        ax.set_title(f"P2 B2-X: |E| — 0 to {vmax:g} MV/cm")
        fig.tight_layout(); fig.savefig(name, dpi=160); plt.close(fig)
        print(name)

    # trench zoom with the neutral boundary
    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    keep = (x[tv].mean(1) >= 2.26) & (x[tv].mean(1) <= 2.48) & \
           (np.abs(y[tv].mean(1)) <= 0.18)
    tp = ax.tripcolor(x, y, tv[keep], facecolors=Emv[keep], shading="flat",
                      cmap="turbo", vmin=0, vmax=5.5)
    ax.tricontour(x, y, tv, u, levels=[1e-8], colors="w", linewidths=1.0)
    for s in (1, -1):
        ax.plot(XS, s*AOF, color="w", lw=1.2, ls="--")
    ax.set_xlim(2.26, 2.48); ax.set_ylim(-0.18, 0.18)
    ax.set_aspect("equal")
    ax.set_xlabel("x (um)"); ax.set_ylabel("y (um)")
    ax.set_title("trench detail: |E| (centroid values), white = neutral boundary")
    fig.colorbar(tp, ax=ax, label="|E| (MV/cm)")
    fig.tight_layout(); fig.savefig("figE3_trench.png", dpi=170); plt.close(fig)
    print("figE3_trench.png")


if __name__ == "__main__":
    main()
