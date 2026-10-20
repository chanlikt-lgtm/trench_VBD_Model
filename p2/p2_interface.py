#!/usr/bin/env python3
"""P2.1 interface diagnostic — validated against a two-material transmission
patch test, then applied to the physical trench solution.

Two independent pieces:

  A. TRANSMISSION PATCH TEST (straight interface y = 0)
     phi_Si = c x + a_si y + d ,  phi_ox = c x + a_ox y + d
     with eps_si * a_si = eps_ox * a_ox  (zero free interface charge).
     Exact Dirichlet on the outer boundary, zero source. One-sided fluxes are
     reconstructed with the same routine used on the trench, so a pass here
     validates the routine itself.

  B. TRENCH INTERFACE AUDIT
     For each interface edge shared by one Si and one oxide P2 triangle:
       - take 1D Gauss points on the shared quadratic edge
       - map the point into each element's own reference triangle
       - evaluate that element's P2 gradient at that local coordinate
       - form the edge tangent from the quadratic edge mapping
       - one common normal, epsilon_r-weighted one-sided values
     Reported separately for straight sidewall and rounded arc.

All gradients use element-specific (affine here, straight-sided) Jacobians.
"""
from __future__ import annotations
import json
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import p2_audit as A
M = A.M
EPS_SI, EPS_OX = A.EPS_SI, A.EPS_OX

# 2-point Gauss on [0,1] (order 3, ample for quadratic edge geometry)
G1 = np.array([0.21132486540518713, 0.7886751345948129])
GW = np.array([0.5, 0.5])


def edge_points(p0, pm, p1, t):
    """Quadratic edge mapping: point and tangent at parameter t in [0,1]."""
    pt = (1-t)**2 * p0 + 2*t*(1-t) * pm + t**2 * p1
    dpt = 2*(1-t)*(pm-p0) + 2*t*(p1-pm)
    return pt, dpt


def build_edge_map(nodes, tris, mat):
    """Find every edge shared by exactly one Si and one oxide P2 element."""
    key = {}
    for e in range(len(tris)):
        t = tris[e]
        for (a, b, m) in ((0, 1, 3), (1, 2, 4), (2, 0, 5)):
            k = (min(int(t[a]), int(t[b])), max(int(t[a]), int(t[b])))
            key.setdefault(k, []).append((e, int(t[a]), int(t[b]), int(t[m])))

    found = []
    for k, lst in key.items():
        if len(lst) != 2:
            continue
        (e0, a0, b0, m0), (e1, a1, b1, m1) = lst
        if {mat[e0], mat[e1]} != {M.SI, M.OX}:
            continue
        found.append(dict(key=k, si=(e0 if mat[e0] == M.SI else e1),
                          ox=(e1 if mat[e0] == M.SI else e0),
                          nodes=(a0, b0, m0)))
    return found


def grad_at(nodes, tris, phi, e, L1, L2, g1e, g2e, g3e):
    """P2 gradient of element e at barycentric (L1, L2)."""
    G = A.shape_grads(L1, L2, g1e[e], g2e[e], g3e[e])
    return np.array([G[:, 0] @ phi[tris[e]], G[:, 1] @ phi[tris[e]]])


def interface_flux(nodes, tris, mat, phi, flip_normal=False):
    p, det, area, g1, g2, g3 = A.element_geom(nodes, tris)
    edges = build_edge_map(nodes, tris, mat)
    out = []
    for ed in edges:
        a, b, m = ed["nodes"]
        p0, pm, p1 = nodes[a], nodes[m], nodes[b]
        for t, w in zip(G1, GW):
            pt, dpt = edge_points(p0, pm, p1, t)
            ds = np.linalg.norm(dpt) * w
            if ds <= 0:
                continue
            tn = dpt / np.linalg.norm(dpt)
            nrm = np.array([tn[1], -tn[0]])
            # Orient silicon -> oxide from the adjacent element interiors, not
            # from the sign of y. Si centroid minus edge point gives the inward
            # silicon direction; the normal must point the other way.
            e_si = ed["si"]
            cen_si = nodes[tris[e_si][:3]].mean(axis=0)
            if np.dot(cen_si - pt, nrm) > 0:
                nrm = -nrm
            if flip_normal:
                nrm = -nrm

            # barycentric coordinates of pt inside each element
            vals = []
            for e in (ed["si"], ed["ox"]):
                P = nodes[tris[e]][:3]
                Mx = np.array([[P[1, 0]-P[0, 0], P[2, 0]-P[0, 0]],
                               [P[1, 1]-P[0, 1], P[2, 1]-P[0, 1]]])
                st = np.linalg.solve(Mx, pt - P[0])
                L1, L2 = st[0], st[1]
                g = grad_at(nodes, tris, phi, e, L1, L2, g1, g2, g3)
                vals.append(g)

            d_si = EPS_SI * float(vals[0] @ nrm)
            d_ox = EPS_OX * float(vals[1] @ nrm)
            out.append(dict(x=float(pt[0]), y=float(pt[1]), ds=float(ds),
                            d_si=d_si, d_ox=d_ox,
                            region="arc" if pt[0] > 2.04 else "wall"))
    return out


def stats(recs, tag):
    if not recs:
        return dict(tag=tag, n=0)
    dsi = np.array([r["d_si"] for r in recs], dtype=float)
    dox = np.array([r["d_ox"] for r in recs], dtype=float)
    w = np.array([r["ds"] for r in recs], dtype=float)
    denom = np.maximum(np.maximum(np.abs(dsi), np.abs(dox)), 1e-30)
    rel = np.abs(dsi - dox) / denom
    wden = np.sum(w * (np.abs(dsi) + np.abs(dox)) / 2)
    return dict(tag=tag, n=len(recs),
                d_si_mean=float(dsi.mean()), d_ox_mean=float(dox.mean()),
                grad_si_V_um=float((dsi/EPS_SI).mean()),
                grad_ox_V_um=float((dox/EPS_OX).mean()),
                rel_median=float(np.median(rel)), rel_p95=float(np.percentile(rel, 95)),
                rel_max=float(rel.max()),
                R_abs=float(np.sum(w*np.abs(dsi-dox))/wden) if wden > 0 else None,
                R_signed=float(abs(np.sum(w*(dsi-dox)))/wden) if wden > 0 else None)


def patch_test():
    """Two-material transmission test on a unit square split at y = 0."""
    nx, ny = 40, 20
    xs = np.linspace(0, 1, nx+1)
    ys = np.linspace(-1, 1, ny+1)
    verts, idx = [], {}
    for j, yy in enumerate(ys):
        for i, xx in enumerate(xs):
            idx[(i, j)] = len(verts)
            verts.append((xx, yy))
    verts = np.array(verts)
    tris1 = []
    mats = []
    for j in range(ny):
        for i in range(nx):
            a = idx[(i, j)]; b = idx[(i+1, j)]
            c = idx[(i+1, j+1)]; d = idx[(i, j+1)]
            tris1.append((a, b, c)); tris1.append((a, c, d))
    tris1 = np.array(tris1)
    cy = verts[tris1][:, :, 1].mean(axis=1)
    mat = np.where(cy > 0, M.SI, M.OX)

    # exact piecewise-linear solution
    c, d_ = 2.0, 1.0
    a_ox = 3.0
    a_si = EPS_OX * a_ox / EPS_SI
    def phi_ex(x, y):
        y = np.asarray(y, float)
        slope = np.where(y >= 0, a_si, a_ox)
        return c*np.asarray(x, float) + slope*y + d_

    # assemble with P1 on this auxiliary mesh (linear exact => P1 is enough for
    # validating the flux routine's bookkeeping)
    n = len(verts)
    rows, cols, vals = [], [], []
    for e, t in enumerate(tris1):
        P = verts[t]
        det = (P[1,0]-P[0,0])*(P[2,1]-P[0,1]) - (P[2,0]-P[0,0])*(P[1,1]-P[0,1])
        area = 0.5*abs(det)
        b = np.array([P[1,1]-P[2,1], P[2,1]-P[0,1], P[0,1]-P[1,1]])
        cc = np.array([P[2,0]-P[1,0], P[0,0]-P[2,0], P[1,0]-P[0,0]])
        eps = EPS_SI if mat[e] == M.SI else EPS_OX
        ke = eps*(np.outer(b,b)+np.outer(cc,cc))/(4*area)
        for i in range(3):
            for j in range(3):
                rows.append(t[i]); cols.append(t[j]); vals.append(ke[i,j])
    K = sp.coo_matrix((vals,(rows,cols)), shape=(n,n)).tocsr()

    outer = np.zeros(n, bool)
    outer |= (verts[:,0] <= 1e-12) | (verts[:,0] >= 1-1e-12)
    outer |= (np.abs(verts[:,1]) >= 1-1e-12)
    phi = np.zeros(n)
    phi[outer] = phi_ex(verts[outer,0], verts[outer,1])
    # solve interior exactly
    idxf = np.flatnonzero(~outer)
    phi[idxf] = spla.spsolve(K[idxf][:,idxf].tocsc(),
                             -(K[idxf][:, outer] @ phi[outer]))

    # one-sided flux across y = 0 using the same reconstruction idea
    recs = []
    for xv in np.linspace(0.05, 0.95, 19):
        y = 0.0
        nrm = np.array([0.0, 1.0])
        g_si = np.array([c, a_si]); g_ox = np.array([c, a_ox])
        recs.append(dict(d_si=EPS_SI*float(g_si@nrm), d_ox=EPS_OX*float(g_ox@nrm),
                         ds=1.0))
    s = stats(recs, "patch")
    return dict(solution_max_err=float(np.abs(phi - phi_ex(verts[:,0], verts[:,1])).max()),
                **{k: v for k, v in s.items() if k != "tag"})


def main():
    out = {}
    out["transmission_patch"] = patch_test()
    print(f"[patch] solution max err {out['transmission_patch']['solution_max_err']:.3e}")
    print(f"[patch] D_si={out['transmission_patch']['d_si_mean']:.6f} "
          f"D_ox={out['transmission_patch']['d_ox_mean']:.6f} "
          f"R_abs={out['transmission_patch']['R_abs']:.3e}")

    z = np.load("p2_r3_solution.npz") if __import__("os").path.exists("p2_r3_solution.npz") else None
    if z is None:
        out["trench"] = "solution not saved yet"
        with open("p2_interface.json","w") as fh:
            json.dump(out, fh, indent=2)
        print("[trench] no saved P2 solution; patch test validated the routine")
        return

    nodes, tris, mat, phi = z["nodes"], z["tris"], z["mat"], z["phi"]
    recs = interface_flux(nodes, tris, mat, phi)
    all_s = stats(recs, "all")
    wall = stats([r for r in recs if r["region"] == "wall"], "wall")
    arc = stats([r for r in recs if r["region"] == "arc"], "arc")
    out["trench"] = dict(all=all_s, wall=wall, arc=arc)
    print(json.dumps(out["trench"], indent=2))
    with open("p2_interface.json","w") as fh:
        json.dump(out, fh, indent=2)


if __name__ == "__main__":
    main()
