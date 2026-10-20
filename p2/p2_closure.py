#!/usr/bin/env python3
"""P2.1 closure checks.

A. Reconcile the r3 wall baseline 0.013879 (all wall) vs 0.011079 (interface
   class): exact numerator/denominator decomposition, not a re-definition.
B. Classify the 2 (r3) + 5 (r4) contact-edge transition samples individually.
C. Independent geometric checks at the rounded tip:
   C1 are the P2 elements straight-sided (midside node at chord midpoint)?
   C2 do the interface vertices lie on the analytic circle?
   C3 chord tangent vs analytic tangent at the Gauss points (angle error)
   C4 Jacobian: det positivity and the 3-point rule's area exactness
   C5 linear-field gradient reproduction through shape_grads + affine Jacobian
   C6 circular-segment area omitted by the polygonal representation
"""
from __future__ import annotations
import json
import numpy as np
import p2_interface as PI
import p2_audit as A

EPS_SI, EPS_OX = A.EPS_SI, A.EPS_OX
R = 0.10
XC = 2.05
A_ARC_END = 0.525


def collect(fn):
    z = np.load(fn)
    nodes, tris, mat, phi, fixed = z["nodes"], z["tris"], z["mat"], z["phi"], z["fixed"]
    p, det, area, g1, g2, g3 = A.element_geom(nodes, tris)
    rows = []
    for ed in PI.build_edge_map(nodes, tris, mat):
        a, b, m = ed["nodes"]
        p0, pm, p1 = nodes[a], nodes[m], nodes[b]
        si_el, ox_el = ed["si"], ed["ox"]
        npin = int(fixed[tris[si_el]].sum())
        allpin = bool(np.all(fixed[tris[si_el]]))
        for t, w in zip(PI.G1, PI.GW):
            pt, dpt = PI.edge_points(p0, pm, p1, t)
            ds = np.linalg.norm(dpt) * w
            tn = dpt / np.linalg.norm(dpt)
            nrm = np.array([tn[1], -tn[0]])
            cen = nodes[tris[si_el][:3]].mean(axis=0)
            if np.dot(cen - pt, nrm) > 0:
                nrm = -nrm
            vals = []
            for e in (si_el, ox_el):
                P = nodes[tris[e]][:3]
                Mx = np.array([[P[1, 0]-P[0, 0], P[2, 0]-P[0, 0]],
                               [P[1, 1]-P[0, 1], P[2, 1]-P[0, 1]]])
                st = np.linalg.solve(Mx, pt - P[0])
                G = A.shape_grads(st[0], st[1], g1[e], g2[e], g3[e])
                vals.append(np.array([G[:, 0] @ phi[tris[e]], G[:, 1] @ phi[tris[e]]]))
            rows.append(dict(x=float(pt[0]), y=float(pt[1]), ds=float(ds),
                             d_si=EPS_SI*float(vals[0] @ nrm),
                             d_ox=EPS_OX*float(vals[1] @ nrm),
                             region="arc" if pt[0] > 2.04 else "wall",
                             si_el=int(si_el), ox_el=int(ox_el), npin=npin,
                             allpin=allpin, edge=(int(a), int(b), int(m))))
    return rows, nodes, tris, mat, phi, fixed


def RD(rows, mask):
    w = np.array([r["ds"] for r, m in zip(rows, mask) if m])
    ds = np.array([r["d_si"] for r, m in zip(rows, mask) if m])
    do = np.array([r["d_ox"] for r, m in zip(rows, mask) if m])
    return (float(np.sum(w*np.abs(ds-do))),
            float(np.sum(w*(np.abs(ds)+np.abs(do))/2)))


def main():
    out = {}
    dref = json.load(open("p2_nearzero.json"))["reference"]["d_ref"]

    # ---------- A + B ----------
    for tag, fn in (("r3", "p2_r3_solution.npz"), ("r4", "p2_r4_solution.npz")):
        rows, nodes, tris, mat, phi, fixed = collect(fn)
        wall = np.array([r["region"] == "wall" for r in rows])
        allpin = np.array([r["allpin"] for r in rows])
        partpin = np.array([(r["npin"] > 0) and not r["allpin"] for r in rows])
        dscale = np.array([(abs(r["d_si"])+abs(r["d_ox"]))/2 for r in rows])
        rel = np.array([abs(r["d_si"]-r["d_ox"])/max(abs(r["d_si"])+abs(r["d_ox"]), 1e-30)
                        for r in rows])
        keep = dscale >= 1e-3*dref

        N_t, D_t = RD(rows, wall)
        N_i, D_i = RD(rows, wall & ~allpin)
        N_e, D_e = RD(rows, wall & allpin)
        print(f"\n[{tag}] wall decomposition")
        print(f"  all        N={N_t:.6f} D={D_t:.6f} R={N_t/D_t:.6f} n={int(wall.sum())}")
        print(f"  interface  N={N_i:.6f} D={D_i:.6f} R={N_i/D_i:.6f} n={int((wall&~allpin).sum())}")
        print(f"  electrode  N={N_e:.6f} D={D_e:.6f} R={N_e/D_e:.6f} n={int((wall&allpin).sum())}")
        print(f"  check (Ni+Ne)/(Di+De) = {(N_i+N_e)/(D_i+D_e):.6f}")
        out[tag] = dict(N_all=N_t, D_all=D_t, N_iface=N_i, D_iface=D_i,
                        N_elec=N_e, D_elec=D_e)

        # transition samples: in interface class, high relative mismatch
        bad = np.flatnonzero(wall & ~allpin & keep & (rel > 0.99))
        print(f"  transition samples (rel>0.99, not fully pinned): {len(bad)}")
        for k in bad:
            r = rows[k]
            el = tris[r["si_el"]]
            pinned = fixed[el]
            print(f"    ({r['x']:.4f},{r['y']:.4f}) npin={r['npin']}/6 "
                  f"d_si={r['d_si']:.4e} d_ox={r['d_ox']:.4e} rel={rel[k]:.3f} "
                  f"phi={np.array2string(phi[el], precision=2)}")
            print(f"       pinned nodes={np.array2string(el[pinned])} "
                  f"dx_to_contact={abs(r['x']-0.54):.4f}")

        # definitive stats with transition samples separated
        trans = wall & partpin
        N_f, D_f = RD(rows, wall & ~allpin & ~trans)
        print(f"  no-free-charge (excl. all-pinned and partial): R={N_f/D_f:.6f} "
              f"n={int((wall&~allpin&~trans).sum())}")
        N_a, D_a = RD(rows, ~wall & ~allpin)
        print(f"  arc no-free-charge: R={N_a/D_a:.6f} n={int((~wall&~allpin).sum())}")
        out[tag]["R_wall_clean"] = N_f/D_f
        out[tag]["n_wall_clean"] = int((wall & ~allpin & ~trans).sum())
        out[tag]["R_arc"] = N_a/D_a

        if tag == "r4":
            # ---------- C. geometry ----------
            print("\n[C] tip geometry")
            p = nodes[tris]
            dev = max(np.abs(p[:, m] - 0.5*(p[:, i] + p[:, j])).max()
                      for (i, j, m) in ((0, 1, 3), (1, 2, 4), (2, 0, 5)))
            print(f"  C1 max |midside - chord midpoint| = {dev:.3e}  "
                  f"({'straight-sided' if dev < 1e-12 else 'curved'})")
            out["C1_mid_dev"] = float(dev)

            det = A.element_geom(nodes, tris)[1]
            print(f"  C4 det>0: {int((det>0).sum())}/{len(det)}  "
                  f"min={det.min():.3e}  max={det.max():.3e}")
            out["C4_min_det"] = float(det.min())

            # interface nodes on the arc
            iface_nodes = sorted({n for ed in PI.build_edge_map(nodes, tris, mat)
                                  for n in ed["nodes"]})
            inp = nodes[iface_nodes]
            sel = (inp[:, 0] >= XC-1e-9) & (inp[:, 0] <= XC+R+1e-9)
            arc_pts = inp[sel]
            rr = np.sqrt((arc_pts[:, 0]-XC)**2 + (np.abs(arc_pts[:, 1])-A_ARC_END)**2)
            print(f"  C2 arc interface nodes={len(arc_pts)}  "
                  f"|r-R| max={np.abs(rr-R).max():.3e}  mean={np.abs(rr-R).mean():.3e}")
            out["C2_radius_dev"] = float(np.abs(rr-R).max())

            # C3 tangent error at Gauss points on arc edges
            ang = []
            seg = 0.0
            nedges = 0
            for ed in PI.build_edge_map(nodes, tris, mat):
                a, b, m = ed["nodes"]
                p0, pm, p1 = nodes[a], nodes[m], nodes[b]
                if not (XC < 0.5*(p0[0]+p1[0]) < XC+R):
                    continue
                nedges += 1
                c = np.linalg.norm(p1-p0)
                th = 2*np.arcsin(min(c/(2*R), 1.0))
                seg += 0.5*R*R*(th - np.sin(th))
                for t in PI.G1:
                    pt, dpt = PI.edge_points(p0, pm, p1, t)
                    tn = dpt/np.linalg.norm(dpt)
                    s = np.sign(pt[1])
                    rad = np.array([pt[0]-XC, abs(pt[1])-A_ARC_END])
                    tan = np.array([-rad[1], rad[0]]); tan /= np.linalg.norm(tan)
                    ang.append(np.degrees(np.arccos(np.clip(abs(np.dot(tn, tan)), -1, 1))))
            print(f"  C3 arc edges={nedges}  tangent angle err: "
                  f"median={np.median(ang):.3f} deg  max={np.max(ang):.3f} deg")
            print(f"  C6 omitted circular segment area total={seg:.6e} um^2 "
                  f"({seg/nedges:.3e} per edge)")
            out["C3_tangent_deg_median"] = float(np.median(ang))
            out["C3_tangent_deg_max"] = float(np.max(ang))
            out["C6_segment_area"] = float(seg)

            # C5 linear-field gradient reproduction on arc-adjacent elements
            arc_els = set()
            for ed in PI.build_edge_map(nodes, tris, mat):
                a, b, m = ed["nodes"]
                p0, p1 = nodes[a], nodes[b]
                if XC < 0.5*(p0[0]+p1[0]) < XC+R:
                    arc_els.add(ed["si"]); arc_els.add(ed["ox"])
            _, _, _, Gg1, Gg2, Gg3 = A.element_geom(nodes, tris)
            errs = []
            for e in arc_els:
                t = tris[e]
                G1 = A.shape_grads(1/3, 1/3, Gg1[e], Gg2[e], Gg3[e])
                gx = np.array([G1[:, 0] @ nodes[t, 0], G1[:, 1] @ nodes[t, 0]])
                gy = np.array([G1[:, 0] @ nodes[t, 1], G1[:, 1] @ nodes[t, 1]])
                errs.append(np.linalg.norm(gx - np.array([1.0, 0.0])))
                errs.append(np.linalg.norm(gy - np.array([0.0, 1.0])))
            print(f"  C5 linear-gradient reproduction on {len(arc_els)} arc elements: "
                  f"max err={np.max(errs):.3e}")
            out["C5_grad_err"] = float(np.max(errs))

    with open("p2_closure.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print("\nwrote p2_closure.json")


if __name__ == "__main__":
    main()
