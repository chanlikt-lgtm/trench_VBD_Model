import json, numpy as np
import p2_front_conditioning as F
from p2_exact_metrics import cut_intervals
from p2_topology_clean import clean
z=np.load("p2_M2_scaled.npz"); nodes,tris,mat,u=z["nodes"],z["tris"],z["mat"],z["u"]
o,_=F.detect_order(nodes,tris)
for x0 in (2.35,2.41):
    ys,us,ds=F.profile(nodes,tris,u,x0,o)
    m,_s=clean(cut_intervals(nodes,tris,mat,u,x0,1e-8,True))
    print(f"x={x0}: intervals {[(round(a,4),round(b,4)) for a,b in m]}")
    for a,b in m:
        for y in (a,b):
            uu=float(np.interp(y,ys,us)); g=float(np.interp(y,ys,ds))
            amb=1e-8/abs(g) if abs(g)>1e-12 else float("inf")
            print(f"   y={y:+.4f}  u={uu:+.3e}  du/dy={g:+.3e}  dy/1e-8={amb:.2e}")
