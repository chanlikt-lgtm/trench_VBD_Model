#!/usr/bin/env python3
"""CP-L2.1/L2.2/L2.3 diagnostics: geometry identity, realized refinement, quality.

Compares S0-L2 against the canonical L0 polygon and against L1 for scaling.
"""
from __future__ import annotations
import hashlib
import json
import numpy as np
from scipy.spatial import cKDTree

import p2_audit as A

M = A.M
LV = "../p2_3_spatial_S0/levels"
BAND_SHA = "D2A00BECA0D38F0183669790CBABC8864D6D467C92C0B609AC85BE0ED62753B3"


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest().upper()


def load(tag):
    z = np.load(f"{LV}/S0_{tag}_p2mesh.npz")
    return z["nodes"], z["tris"], z["mat"]


def boundary_edges(nodes, tris, mat):
    sel = tris[mat == M.SI][:, :3]
    e = np.sort(np.vstack([sel[:, [0, 1]], sel[:, [1, 2]], sel[:, [2, 0]]]), axis=1)
    uniq, cnt = np.unique(e, axis=0, return_counts=True)
    b = uniq[cnt == 1]
    return nodes[b[:, 0]], nodes[b[:, 1]]


def pt_seg_dist(pts, Aseg, Bseg):
    best = np.full(len(pts), np.inf)
    for i in range(len(Aseg)):
        a = Aseg[i]; ab = Bseg[i] - a; L2 = float(ab @ ab)
        t = np.clip(((pts - a) @ ab) / L2, 0.0, 1.0)
        d = np.linalg.norm(pts - (a + t[:, None] * ab), axis=1)
        best = np.minimum(best, d)
    return best


def midside_dev(nodes, tris):
    P = nodes[tris]
    return float(max(np.abs(P[:, 3] - 0.5*(P[:, 0]+P[:, 1])).max(),
                     np.abs(P[:, 4] - 0.5*(P[:, 1]+P[:, 2])).max(),
                     np.abs(P[:, 5] - 0.5*(P[:, 2]+P[:, 0])).max()))


def quality(nodes, tris, mat):
    P = nodes[tris][:, :3]
    det = ((P[:, 1, 0]-P[:, 0, 0])*(P[:, 2, 1]-P[:, 0, 1]) -
           (P[:, 2, 0]-P[:, 0, 0])*(P[:, 1, 1]-P[:, 0, 1]))
    area = 0.5*np.abs(det)
    edges = np.stack([np.linalg.norm(P[:, 1]-P[:, 0], axis=1),
                      np.linalg.norm(P[:, 2]-P[:, 1], axis=1),
                      np.linalg.norm(P[:, 0]-P[:, 2], axis=1)], 1)
    emax = edges.max(1); hmin = 2*area/np.maximum(emax, 1e-30)
    aspect = emax/np.maximum(hmin, 1e-30)
    ang = []
    for k in range(3):
        v1 = P[:, (k+1) % 3]-P[:, k]; v2 = P[:, (k+2) % 3]-P[:, k]
        cs = np.sum(v1*v2, 1)/(np.linalg.norm(v1, axis=1)*np.linalg.norm(v2, axis=1)+1e-30)
        ang.append(np.degrees(np.arccos(np.clip(cs, -1, 1))))
    ang = np.concatenate(ang)
    orient = {}
    for code in np.unique(mat):
        m = mat == code
        orient[str(int(code))] = dict(n=int(m.sum()), pos=int((det[m] > 0).sum()),
                                      neg=int((det[m] < 0).sum()),
                                      zero=int((det[m] == 0).sum()))
    return dict(min_angle=float(ang.min()), p01=float(np.percentile(ang, 1)),
                aspect_max=float(aspect.max()), aspect_p99=float(np.percentile(aspect, 99)),
                n_gt3=int((aspect > 3).sum()), n_gt5=int((aspect > 5).sum()),
                zero_area=int((area <= 0).sum()), min_area=float(area.min()),
                orientation=orient)


def regions(nodes, tris, mat, band):
    P = nodes[tris][:, :3]
    area = 0.5*np.abs(((P[:, 1, 0]-P[:, 0, 0])*(P[:, 2, 1]-P[:, 0, 1]) -
                       (P[:, 2, 0]-P[:, 0, 0])*(P[:, 1, 1]-P[:, 0, 1])))
    cx = P[:, :, 0].mean(1); cy = P[:, :, 1].mean(1)
    si = mat == M.SI
    h = np.sqrt(4*area/np.sqrt(3))
    d, _ = cKDTree(band).query(np.column_stack([cx, cy]))
    out = {}
    # distance to the two canonical arc circles (radius R about the arc centres)
    R = abs(M.A_BOT - M.A_ARC_END)
    darc = np.minimum(np.abs(np.hypot(cx - M.TIP_START, cy - M.A_ARC_END) - R),
                      np.abs(np.hypot(cx - M.TIP_START, cy + M.A_ARC_END) - R))
    defs = {"arc<=0.01": darc <= 0.01,
            "band<=0.04": d <= 0.04,
            "trans0.04-0.06": (d > 0.04) & (d <= 0.06),
            "tip": (cx > 2.30) & (cx < 2.42) & (np.abs(cy) < 0.02),
            "bodycut": np.abs(cx - A.BODY_BOTTOM) <= 0.05,
            "background": d > 0.15}
    for k, m in defs.items():
        m = m & si
        if m.sum() > 5:
            out[k] = dict(n=int(m.sum()), p50=float(np.median(h[m])),
                          p90=float(np.percentile(h[m], 90)))
    return out


def main():
    band = np.load("../p2_3_spatial_S0/levels/../p2_band_contour.npz")["pts"] \
        if False else np.load("p2_band_contour.npz")["pts"]
    out = {"band_sha": sha("p2_band_contour.npz"), "band_sha_expected": BAND_SHA,
           "band_sha_match": sha("p2_band_contour.npz") == BAND_SHA}
    data = {t: load(t) for t in ("L0", "L1", "L2")}
    out["mesh_sha"] = {t: sha(f"{LV}/S0_{t}_p2mesh.npz") for t in data}

    for t in ("L1", "L2"):
        n, tr, m = data[t]
        P = n[tr][:, :3]
        area = 0.5*np.abs(((P[:, 1, 0]-P[:, 0, 0])*(P[:, 2, 1]-P[:, 0, 1]) -
                           (P[:, 2, 0]-P[:, 0, 0])*(P[:, 1, 1]-P[:, 0, 1])))
        si = m == M.SI
        print(f"[{t}] nodes {len(n)} elems {len(tr)} Si {area[si].sum():.12f} "
              f"midside {midside_dev(n, tr):.2e}", flush=True)
        out[t] = dict(nodes=int(len(n)), elements=int(len(tr)),
                      si_area=float(area[si].sum()),
                      midside_dev=midside_dev(n, tr),
                      quality=quality(n, tr, m),
                      regions=regions(n, tr, m, band))

    # geometry identity of L2 boundary against the L0 polygon
    n0, t0, m0 = data["L0"]
    n2, t2, m2 = data["L2"]
    Aseg, Bseg = boundary_edges(n0, t0, m0)
    s2 = t2[m2 == M.SI][:, :3]
    e2 = np.sort(np.vstack([s2[:, [0, 1]], s2[:, [1, 2]], s2[:, [2, 0]]]), axis=1)
    u2, c2 = np.unique(e2, axis=0, return_counts=True)
    pts = n2[np.unique(u2[c2 == 1].ravel())]
    d = pt_seg_dist(pts, Aseg, Bseg)
    out["L2_boundary_drift"] = dict(n=int(len(pts)), max=float(d.max()),
                                    p99=float(np.percentile(d, 99)))
    print(f"[L2 vs L0 polygon] boundary nodes {len(pts)}  max {d.max():.3e} um  "
          f"p99 {np.percentile(d,99):.3e}", flush=True)

    # scaling ratios L2/L1
    rat = {}
    for k in out["L2"]["regions"]:
        if k in out["L1"]["regions"]:
            a = out["L1"]["regions"][k]; b = out["L2"]["regions"][k]
            rat[k] = dict(p50=round(b["p50"]/a["p50"], 4), p90=round(b["p90"]/a["p90"], 4))
    out["ratio_L2_over_L1"] = rat
    out["target_ratio"] = 2**-0.5
    print("[ratios L2/L1] target 0.7071:", json.dumps(rat), flush=True)
    json.dump(out, open("../p2_3_spatial_S0/CP-L2_geometry.json", "w"), indent=2)
    print("wrote CP-L2_geometry.json")


if __name__ == "__main__":
    main()
