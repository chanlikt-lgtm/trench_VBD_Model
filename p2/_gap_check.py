import numpy as np
import p2_front_conditioning as F
from p2_exact_metrics import cut_intervals
from p2_topology_clean import clean

MESH=[("r3","p2_p3_r3_polished.npz"),("r4","p2_p3_r4_polished.npz"),
      ("mid","p2_p3_rmid_polished.npz"),("L1b","p2_L1b_scaled.npz"),("L1","p2_L1_scaled.npz")]

print("cut  mesh  eps      n_raw  intervals                                        gap_max(V)      n_clean")
for x0 in (2.35,2.41):
    for n,f in MESH:
        z=np.load(f); o,_=F.detect_order(z["nodes"],z["tris"])
        ys,us,ds=F.profile(z["nodes"],z["tris"],z["u"],x0,o)
        for eps in (1e-8,1e-6,1e-5):
            iv=cut_intervals(z["nodes"],z["tris"],z["mat"],z["u"],x0,eps,True)
            iv=sorted([(float(a),float(b)) for a,b in iv])
            gaps=[]
            for (a1,b1),(a2,b2) in zip(iv[:-1],iv[1:]):
                m=(ys>b1)&(ys<a2)
                gaps.append(float(us[m].max()) if m.any() else float("nan"))
            mc,_sc=clean(iv)
            gtxt=",".join(f"{g:.2e}" for g in gaps) if gaps else "-"
            itxt=" ".join(f"({a:+.4f},{b:+.4f})" for a,b in iv)
            print(f"{x0:.2f} {n:4s} {eps:g}  {len(iv):5d}  {itxt:46s} {gtxt:14s} {len(mc):5d}")
    print()
