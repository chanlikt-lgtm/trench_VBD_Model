import numpy as np
import p2_audit as A
from p2_bernstein import B
M=A.M; EPS=1e-8
def bands(fn):
    z=np.load(fn); nodes,tris,mat,u=z["nodes"],z["tris"],z["mat"],z["u"]
    area=np.abs(A.element_geom(nodes,tris)[2]); P=nodes[tris][:,:3]
    cx=P[:,:,0].mean(1); cy=P[:,:,1].mean(1)
    ch=(mat==M.SI)&(cx>A.BODY_BOTTOM+1e-12)
    c=(B@u[tris].T).T-EPS
    mixed=ch&~((c<=0).all(1))&~((c>0).all(1))
    h=np.sqrt(4*area/np.sqrt(3))
    out={}
    for a in np.arange(2.20,2.50,0.05):
        m=mixed&(cx>=a)&(cx<a+0.05)
        if m.sum()>=5: out[round(a,2)]=(int(m.sum()),float(np.median(h[m])))
    return out
L1b=bands("p2_L1b_scaled.npz"); L1=bands("p2_L1_scaled.npz"); M2=bands("p2_M2_scaled.npz")
print("x band   L1b(n,h)        M2(n,h)         L1(n,h)")
for k in sorted(set(L1b)|set(L1)):
    f=lambda d: f"{d[k][0]:3d},{d[k][1]:.5f}" if k in d else "   -,-     "
    print(f"{k:.2f}   {f(L1b):14s} {f(M2):14s} {f(L1):14s}")
