#!/usr/bin/env python3
"""CP-S0.2 fixture suite: donor quadrature, area estimator, front extractor.

All fixtures are self-contained; the B2 mesh is used only for the global
charge-conservation check.
"""
from __future__ import annotations
import json
import numpy as np

import p2_audit as A
import p2_area_bounds as AB
import p2_exact_donor as D
from p2_front_conditioning import p2_eval

M = A.M
XB = A.BODY_BOTTOM
RES = {}


def tri_mesh(offset_x=0.0, at=0.7):
    """Unit right triangle with area 0.5, placed inside the charged region."""
    nodes = np.array([[at, 0.0], [at+1.0, 0.0], [at, 1.0],
                      [at+0.5, 0.0], [at+0.5, 0.5], [at, 0.5]])
    tris = np.array([[0, 1, 2, 3, 4, 5]])
    mat = np.array([M.SI])
    return nodes, tris, mat


def fixture_donor():
    nodes, tris, mat = tri_mesh()
    area = abs(A.element_geom(nodes, tris)[2][0])

    # 1. fully charged: identity rho*A*[0,0,0,1/3,1/3,1/3]
    Pv = nodes[tris[0][:3]]
    I = D.int_N(Pv, -1e9, area)
    want = area*np.array([0, 0, 0, 1/3, 1/3, 1/3])
    e1 = float(np.abs(I - want).max())

    # 2. fully body-side: zero
    I2 = D.int_N(Pv, 1e9, area)
    e2 = float(np.abs(I2).max())

    # 3. cut element: sum equals clipped area, all components >= 0
    P = Pv
    xb = 0.7 + 0.4                       # cuts the triangle
    I3 = D.int_N(P, xb, area)
    from p2_charge_mask_audit import clip_area
    Aclip = clip_area([tuple(q) for q in P], xb)
    e3 = abs(I3.sum() - Aclip)
    neg = float(I3.min())

    RES["donor"] = dict(full_identity_max_err=e1, body_side_max=I2.max(),
                        cut_sum_err=float(e3), cut_min_component=neg,
                        cut_area=Aclip,
                        note="negative vertex components are expected: N=L(2L-1) is sign-changing",
                        passed=bool(e1 < 1e-12 and I2.max() == 0.0
                                    and e3 < 1e-12 and np.isfinite(I3).all()))


def fixture_area():
    nodes, tris, mat = tri_mesh()
    area = abs(A.element_geom(nodes, tris)[2][0])

    # (a) u = 0 everywhere -> whole triangle neutral
    u0 = np.zeros(len(nodes))
    lo, up, _, _ = AB.bounds(nodes, tris, mat, u0)
    a_ok = abs(lo - area) < 1e-15 and abs(up - area) < 1e-15

    # (b) u = +1 everywhere -> nothing neutral
    u1 = np.ones(len(nodes))
    lo1, up1, _, _ = AB.bounds(nodes, tris, mat, u1)
    b_ok = abs(lo1) < 1e-15 and abs(up1) < 1e-15

    # (c) manufactured linear: u = L1 = x - at, eps = 0.5 -> area 0.375 of 0.5
    #     (documented fixture from the F2 record)
    xs = nodes[:, 0] - nodes[0, 0]
    u_lin = xs
    old_eps = AB.EPS
    AB.EPS = 0.5
    lo2, up2, _, _ = AB.bounds(nodes, tris, mat, u_lin)
    AB.EPS = old_eps
    c_ok = (abs(lo2 - 0.375) < 1e-12) and (0.0 <= up2 - lo2 < 1e-5)

    RES["area"] = dict(zero_field=[lo, up], zero_ok=bool(a_ok),
                       positive_field=[lo1, up1], positive_ok=bool(b_ok),
                       linear_fixture=[float(lo2), float(up2)],
                       linear_ok=bool(c_ok),
                       passed=bool(a_ok and b_ok and c_ok))


def fixture_front():
    P = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    V = np.array([[0., 0.], [1., 0.], [0., 1.], [.5, 0.], [.5, .5], [0., .5]])
    ue = 2*V[:, 0] - 3*V[:, 1]
    order = ((0, 1), (1, 2), (2, 0))
    errs = []
    for (x, y) in ((0.25, 0.25), (0.6, 0.2), (0.1, 0.7)):
        uu, dd, _ = p2_eval(P, ue, x, np.array([y]), order)
        errs.append((abs(float(uu[0]) - (2*x - 3*y)), abs(float(dd[0]) + 3.0)))
    RES["front"] = dict(max_value_err=float(max(e[0] for e in errs)),
                        max_grad_err=float(max(e[1] for e in errs)),
                        passed=bool(max(max(e) for e in errs) < 1e-12))


def fixture_global():
    z = np.load("p2_B2_p2mesh.npz")
    _, tot = D.exact_load(z["nodes"], z["tris"], z["mat"])
    m = json.load(open("p2_B2_mesh.json"))
    RES["global"] = dict(sum_f_over_rho=tot, exact_charged=m["exact_charged"],
                         match=bool(abs(tot - m["exact_charged"]) < 1e-12))


def main():
    fixture_donor()
    fixture_area()
    fixture_front()
    fixture_global()
    ok = (RES["donor"]["passed"] and RES["area"]["passed"]
          and RES["front"]["passed"] and RES["global"]["match"])
    RES["all_passed"] = bool(ok)
    print(json.dumps(RES, indent=2))
    json.dump(RES, open("p2_s0_tests.json", "w"), indent=2)
    print("\nCP-S0.2:", "PASS" if ok else "FAIL")


if __name__ == "__main__":
    main()
