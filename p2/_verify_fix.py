import numpy as np
import p2_front_conditioning as F
from p2_exact_metrics import cut_intervals
from p2_topology_clean import clean

# --- synthetic check through the patched evaluator ---
P=np.array([[0.,0.],[1.,0.],[0.,1.]])
V=np.array([[0.,0.],[1.,0.],[0.,1.],[.5,0.],[.5,.5],[0.,.5]])
ue=2*V[:,0]-3*V[:,1]
o=((0,1),(1,2),(2,0))
for (x,y) in ((0.25,0.25),(0.6,0.2),(0.1,0.7)):
    uu,dd,_=F.p2_eval(P,ue,x,np.array([y]),o)
    print(f"  ({x:.2f},{y:.2f}) exact {2*x-3*y:+.6f}  patched {float(uu[0]):+.6f}  dudy {float(dd[0]):+.4f}")

MESH=[("r3","p2_p3_r3_polished.npz"),("r4","p2_p3_r4_polished.npz"),
      ("mid","p2_p3_rmid_polished.npz"),("L1b","p2_L1b_scaled.npz"),("L1","p2_L1_scaled.npz")]
for x0 in (2.35,2.41):
    print(f"\n=== x={x0} (corrected extractor) ===")
    for n,f in MESH:
        z=np.load(f)
        row=[]
        for eps in (1e-8,1e-6,1e-5):
            m,_s=clean(cut_intervals(z["nodes"],z["tris"],z["mat"],z["u"],x0,eps,True))
            row.append(" ".join(f"({a:+.4f},{b:+.4f})" for a,b in m))
        print(f"  {n:4s} 1e-8 {row[0]}")
        print(f"       1e-6 {row[1]}")
        print(f"       1e-5 {row[2]}")

