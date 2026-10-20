import numpy as np
import p2_front_conditioning as F
import p2_audit as A

M=A.M
def eval_pt(nodes,tris,mat,u,x0,y,order):
    P=nodes[tris]
    xmin=P[:,:,0].min(1); xmax=P[:,:,0].max(1)
    cx=P[:,:,:2].mean(1)
    sel=np.flatnonzero((xmin<=x0+1e-12)&(xmax>=x0-1e-12)&(mat==M.SI))
    best=None
    for e in sel:
        p=P[e,:3]
        out=F.p2_eval(p,u[tris[e]],x0,np.array([y]),order)
        if out is None: continue
        uu,dd,yy=out
        if len(uu)==0: continue
        # inside test: p2_eval returns only inside samples
        v=float(uu[0]); g=float(dd[0])
        if best is None or abs(v-1e-8)<abs(best[0]-1e-8): best=(v,g)
    return best

for name,f in (("L1b","p2_L1b_scaled.npz"),("L1","p2_L1_scaled.npz")):
    z=np.load(f); nodes,tris,mat,u=z["nodes"],z["tris"],z["mat"],z["u"]
    order,_=F.detect_order(nodes,tris)
    print(f"\n=== {name} at x=2.35 ===")
    for y in (-0.10,-0.09,-0.085,-0.080,-0.0774,-0.075,-0.0718,-0.068,-0.060,-0.050):
        r=eval_pt(nodes,tris,mat,u,2.35,y,order)
        if r is None: print(f"  y={y:+.4f}  (no Si element)"); continue
        print(f"  y={y:+.4f}  u={r[0]:+.6e}  du/dy={r[1]:+.4e}")
