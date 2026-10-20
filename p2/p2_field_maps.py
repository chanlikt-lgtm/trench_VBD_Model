#!/usr/bin/env python3
"""Three 2D field maps on the validated intermediate P2 solution.

1. potential phi(x,y) with material boundaries (Si, SiO2, G, SG)
2. |E|(x,y) with streamlines; unclipped and clipped scales
3. neutral-region and depletion-front contours on the silicon geometry
"""
from __future__ import annotations
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree

import p2_audit as A
import mesh_stage2b as M2

M = A.M
V_D = A.V_DRAIN


def main():
    z = np.load("p2_p3_rmid_polished.npz")
    nodes, tris, mat, phi, u = z["nodes"], z["tris"], z["mat"], z["phi"], z["u"]
    tv = tris[:, :3]
    x, y = nodes[:, 0], nodes[:, 1]

    # element centroid field
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
    Eabs = np.linalg.norm(E, axis=1)

    # material boundaries
    xs = np.linspace(0.0, 2.42, 400)
    a = M2.a_of(xs)

    def draw_bounds(ax):
        for s in (1, -1):
            ax.plot(xs, s*a, color="k", lw=0.8, alpha=0.8)
            ax.plot(xs, s*(a + 0.10), color="k", lw=0.5, ls="--", alpha=0.6)
        ax.set_xlim(0, 2.42); ax.set_ylim(-0.675, 0.675)
        ax.set_aspect("equal"); ax.set_xlabel("x (um)"); ax.set_ylabel("y (um)")

    # ---- map 1: potential ----
    fig, ax = plt.subplots(figsize=(10, 5.2))
    tcf = ax.tricontourf(x, y, tv, phi, levels=40, cmap="viridis")
    draw_bounds(ax)
    cb = fig.colorbar(tcf, ax=ax, label="phi (V)")
    ax.text(1.6, 0.0, "Si", ha="center", va="center", color="w", fontsize=10)
    ax.text(1.6, 0.60, "SiO2", ha="center", va="center", color="w", fontsize=9)
    ax.text(0.32, 0.55, "G", ha="center", va="center", color="w", fontsize=9)
    ax.text(1.35, 0.60, "SG", ha="center", va="center", color="w", fontsize=9)
    ax.set_title("P2 intermediate: potential (V_D - u)")
    fig.tight_layout(); fig.savefig("p2_mid_phi.png", dpi=150); plt.close(fig)

    # ---- map 2: |E|, two scales + streamlines ----
    # coarse grid for streamlines, KD-tree locator over element centroids
    cx = nodes[tv][:, :, 0].mean(axis=1)
    cy = nodes[tv][:, :, 1].mean(axis=1)
    tree = cKDTree(np.column_stack([cx, cy]))
    gx = np.linspace(0.01, 2.41, 120)
    gy = np.linspace(-0.66, 0.66, 60)
    GX, GY = np.meshgrid(gx, gy)
    U = np.full(GX.shape, np.nan); V = np.full(GX.shape, np.nan)
    for i in range(GX.shape[0]):
        pts = np.column_stack([GX[i], GY[i]])
        _, idxs = tree.query(pts, k=12)
        for j in range(len(pts)):
            found = False
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
                    found = True
                    break
            if not found:
                U[i, j] = V[i, j] = np.nan
    Um = np.nan_to_num(U); Vm = np.nan_to_num(V)

    fig, axes = plt.subplots(2, 1, figsize=(10, 9.5), sharex=True, sharey=True)
    for ax, vmax, label in ((axes[0], float(Eabs.max()), "unclipped"),
                            (axes[1], 100.0, "clipped at 100 V/um")):
        tp = ax.tripcolor(x, y, tv, facecolors=Eabs, shading="flat",
                          cmap="inferno", vmin=0, vmax=vmax)
        ax.streamplot(GX, GY, Um, Vm, color="w", linewidth=0.4,
                      density=0.8, arrowsize=0.5)
        draw_bounds(ax)
        fig.colorbar(tp, ax=ax, label="|E| (V/um)")
        ax.set_title(f"|E|, {label}")
    fig.tight_layout(); fig.savefig("p2_mid_E.png", dpi=150); plt.close(fig)

    # ---- map 3: neutral and front contours ----
    fig, ax = plt.subplots(figsize=(10, 5.2))
    simask = mat == M.SI
    face = np.where(simask, 0.85, 0.97)
    ax.tripcolor(x, y, tv, facecolors=face, shading="flat",
                 cmap="gray", vmin=0, vmax=1)
    levels = [1e-8, 1e-4, 1e-2, 1e-1, 1.0]
    cs = ax.tricontour(x, y, tv, u, levels=levels, colors="k", linewidths=0.8)
    ax.clabel(cs, fmt=lambda v: f"{v:g}", fontsize=7)
    front = ax.tricontour(x, y, tv, u, levels=[1e-8], colors="r", linewidths=2.0)
    draw_bounds(ax)
    handles = [Line2D([0], [0], color="r", lw=2, label="u = 1e-8 (front)"),
               Line2D([0], [0], color="k", lw=0.8, label="u contours")]
    ax.legend(handles=handles, loc="upper right", fontsize=8)
    ax.set_title("Neutral region and depletion front (grey = Si)")
    fig.tight_layout(); fig.savefig("p2_mid_neutral.png", dpi=150); plt.close(fig)
    print("wrote p2_mid_phi.png, p2_mid_E.png, p2_mid_neutral.png")
    print(f"|E| max at element centroids = {Eabs.max():.4f} V/um "
          f"(not a converged tip maximum)")


if __name__ == "__main__":
    main()
