#!/usr/bin/env python3
"""P2.3 post-processing: comparable neutral-set metrics and front location."""
from __future__ import annotations
import json
import numpy as np
import p2_audit as A

M = A.M
V_D = A.V_DRAIN


def main():
    z = np.load("p2_p3_r3.npz")
    nodes, tris, mat = z["nodes"], z["tris"], z["mat"]
    phi, u, active = z["phi"], z["u"], z["active"]
    area = np.abs(A.element_geom(nodes, tris)[2])
    si = mat == M.SI
    tol = 1e-8

    any_node = si & (np.abs(u[tris]) < tol).any(axis=1)
    all_nodes = si & (np.abs(u[tris]) < tol).all(axis=1)
    act_rows = si & active

    print(f"[neutral] element-level Bernstein active area = {area[act_rows].sum():.9f} "
          f"({int(act_rows.sum())} elements)")
    print(f"[neutral] any node |u|<{tol:g}            = {area[any_node].sum():.9f} "
          f"({int(any_node.sum())} elements)")
    print(f"[neutral] all nodes |u|<{tol:g}           = {area[all_nodes].sum():.9f} "
          f"({int(all_nodes.sum())} elements)")
    print(f"[neutral] charged Si area                = {area[si & ~(nodes[tris][:, :3, 0].mean(axis=1) <= A.BODY_BOTTOM)].sum():.9f}")

    nodes_si = np.unique(tris[si])
    uu = u[nodes_si]
    near = nodes_si[np.abs(uu) < 1e-6]
    print(f"[neutral] Si nodes with |u|<1e-6: {len(near)}")
    if len(near):
        print(f"          x range {nodes[near,0].min():.4f}..{nodes[near,0].max():.4f} "
              f"y range {nodes[near,1].min():.4f}..{nodes[near,1].max():.4f}")

    # front: furthest downstream point with u < 1e-6 (Si, free nodes)
    free_guess = np.abs(u) < 1e-6
    xs = nodes[nodes_si][np.abs(u[nodes_si]) < 1e-6, 0]
    print(f"[front] max x with |u|<1e-6 = {xs.max():.5f}" if len(xs) else "[front] none")

    # phi overshoot check
    print(f"[phi] max = {phi.max():.9f}  (V_D = {V_D:.9f})  overshoot = {phi.max()-V_D:+.3e}")

    out = dict(
        bernstein_active_area=float(area[act_rows].sum()),
        bernstein_active_elements=int(act_rows.sum()),
        any_node_area=float(area[any_node].sum()),
        all_nodes_area=float(area[all_nodes].sum()),
        near_zero_nodes=int(len(near)),
        max_x_near_zero=float(xs.max()) if len(xs) else None,
        phi_max=float(phi.max()),
        overshoot=float(phi.max()-V_D),
    )
    json.dump(out, open("p2_p3_post.json", "w"), indent=2)
    print("wrote p2_p3_post.json")


if __name__ == "__main__":
    main()
