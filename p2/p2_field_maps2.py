#!/usr/bin/env python3
"""Field maps styled to match the reference figure.

|E| in MV/cm with a 0-5.5 scale (reference) plus a 0-2.0 bulk-contrast panel;
thick white dashed Si/SiO2 interfaces; white streamlines; labelled regions.
"""
from __future__ import annotations
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree

import p2_audit as A
import mesh_stage2b as M2

M = A.M


def main():
    z = np.load("p2_p3_rmid_polished.npz")
    nodes, tris, mat, phi, u = z["nodes"], z["tris"], z["mat"], z["phi"], z["u"]
    tv = tris[:, :3]
    x, y = nodes[:, 0], nodes[:, 1]

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
    Emv = np.linalg.norm(E, axis=1)/100.0     # V/um -> MV/cm

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

    # ---- potential ----
    fig, ax = plt.subplots(figsize=(10, 5.2))
    tcf = ax.tricontourf(x, y, tv, phi, levels=40, cmap="viridis")
    bounds(ax, color="w")
    labels(ax)
    fig.colorbar(tcf, ax=ax, label="phi (V)")
    ax.set_title("P2 intermediate: potential")
    fig.tight_layout(); fig.savefig("p2_mid_phi_v2.png", dpi=160); plt.close(fig)

    # ---- |E| ----
    cx = nodes[tv][:, :, 0].mean(axis=1)
    cy = nodes[tv][:, :, 1].mean(axis=1)
    tree = cKDTree(np.column_stack([cx, cy]))
    gx = np.linspace(0.01, 2.41, 130)
    gy = np.linspace(-0.66, 0.66, 66)
    GX, GY = np.meshgrid(gx, gy)
    U = np.full(GX.shape, np.nan); V = np.full(GX.shape, np.nan)
    for i in range(GX.shape[0]):
        pts = np.column_stack([GX[i], GY[i]])
        _, idxs = tree.query(pts, k=12)
        for j in range(len(pts)):
            for e in np.atleast_1d(idxs[j]):
                P = nodes[tv[e]]
                Mx = np.array([[P[1, 0]-P[0, 0], P[2, 0]-P[0, 0]],
                               [P[1, 1]-P[0, 1], P[2, 1]-P[0, 1]]])
                try:
                    st = np.linalg.solve(Mx, pts[j]-P[0])
                except np.linalg.LinAlgError:
                    continue
                l1, l2 = st
                if l1 >= -1e-9 and l2 >= -1e-9 and (1-l1-l2) >= -1e-9:
                    U[i, j], V[i, j] = E[e]
                    break
    Um = np.nan_to_num(U); Vm = np.nan_to_num(V)

    fig, axes = plt.subplots(2, 1, figsize=(10, 9.6), sharex=True, sharey=True)
    for ax, vmax, title in ((axes[0], 5.5, "0-5.5 MV/cm (reference scale)"),
                            (axes[1], 2.0, "0-2.0 MV/cm (bulk contrast)")):
        tp = ax.tripcolor(x, y, tv, facecolors=Emv, shading="flat",
                          cmap="turbo", vmin=0, vmax=vmax)
        ax.streamplot(GX, GY, Um, Vm, color="w", linewidth=0.5,
                      density=0.9, arrowsize=0.55)
        bounds(ax, color="w")
        labels(ax)
        fig.colorbar(tp, ax=ax, label="|E| (MV/cm)")
        ax.set_title(title)
    fig.suptitle("P2 intermediate: |E| with Si/SiO2 interfaces (white dashed)")
    fig.tight_layout(); fig.savefig("p2_mid_E_v2.png", dpi=160); plt.close(fig)

    # ---- neutral / front ----
    fig, ax = plt.subplots(figsize=(10, 5.2))
    simask = mat == M.SI
    face = np.where(simask, 0.82, 0.97)
    ax.tripcolor(x, y, tv, facecolors=face, shading="flat", cmap="gray", vmin=0, vmax=1)
    levels = [1e-8, 1e-4, 1e-2, 1e-1, 1.0]
    cs = ax.tricontour(x, y, tv, u, levels=levels, colors="k", linewidths=0.8)
    ax.clabel(cs, fmt=lambda v: f"{v:g}", fontsize=7)
    ax.tricontour(x, y, tv, u, levels=[1e-8], colors="r", linewidths=2.2)
    bounds(ax, color="k", lw=1.2)
    labels(ax)
    handles = [Line2D([0], [0], color="r", lw=2.2, label="u = 1e-8 (front)"),
               Line2D([0], [0], color="k", lw=0.8, label="u contours")]
    ax.legend(handles=handles, loc="upper right", fontsize=8)
    ax.set_title("Neutral region and depletion front (grey = Si)")
    fig.tight_layout(); fig.savefig("p2_mid_neutral_v2.png", dpi=160); plt.close(fig)

    print("wrote p2_mid_phi_v2.png, p2_mid_E_v2.png, p2_mid_neutral_v2.png")
    print(f"|E| max = {Emv.max():.3f} MV/cm (element centroids; not a converged tip maximum)")


if __name__ == "__main__":
    main()
