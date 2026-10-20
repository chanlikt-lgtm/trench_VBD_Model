import numpy as np, json
import p2_front_conditioning as F
from p2_exact_metrics import cut_intervals
from p2_topology_clean import clean

MESH=[("mid","p2_p3_rmid_polished.npz"),("L1b","p2_L1b_scaled.npz"),("L1","p2_L1_scaled.npz")]
res={}
for x0 in (2.35,2.41):
    prof={}; ends={}
    for n,f in MESH:
        z=np.load(f); o,_=F.detect_order(z["nodes"],z["tris"])
        ys,us,ds=F.profile(z["nodes"],z["tris"],z["u"],x0,o)
        prof[n]=(ys,us,ds)
        iv=cut_intervals(z["nodes"],z["tris"],z["mat"],z["u"],x0,1e-8,True)
        m,_s=clean(iv)
        ends[n]=[float(y) for ab in m for y in ab]
    print(f"--- x={x0}: du at each endpoint, predicted vs observed shift ---")
    print(f"{'A':>4s} {'B':>4s} {'y*':>9s} {'du(y*)':>10s} {'(du/dy)_B':>10s} {'pred dy':>10s} {'obs dy':>9s}")
    rows=[]
    for a,b in (("L1b","L1"),("mid","L1b")):
        ya,ua,da=prof[a]; yb,ub,db=prof[b]
        for y in ends[a]:
            du=float(np.interp(y,yb,ub))-float(np.interp(y,ya,ua))
            g=float(np.interp(y,yb,db))
            pred=du/g if abs(g)>1e-12 else float("inf")
            near=min(ends[b], key=lambda t: abs(t-y))
            obs=near-y
            print(f"{a:>4s} {b:>4s} {y:+9.4f} {du:+10.2e} {g:+10.3e} {pred:+10.4f} {obs:+9.4f}")
            rows.append(dict(A=a,B=b,y=y,du=du,dudy=g,pred=pred,obs=obs))
    res[f"{x0}"]=rows
json.dump(res,open("p2_front_amp.json","w"),indent=2)
print("wrote p2_front_amp.json")
