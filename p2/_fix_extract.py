import numpy as np
import p2_audit as A
from p2_exact_metrics import cut_intervals
from p2_topology_clean import clean
M=A.M

def u_at(P, ue, x0, y):
    Mx=np.array([[P[1,0]-P[0,0],P[2,0]-P[0,0]],[P[1,1]-P[0,1],P[2,1]-P[0,1]]])
    Mi=np.linalg.inv(Mx)
    a,b=Mi@(np.array([x0,y])-P[0])
    L1=1.0-a-b; L2=a; L3=b
    N=np.array([L1*(2*L1-1),L2*(2*L2-1),L3*(2*L3-1),4*L1*L2,4*L2*L3,4*L3*L1])
    return float(N@ue), (L1,L2,L3)

def corrected(nodes,tris,mat,u,x0,eps,Ns=129):
    P=nodes[tris]; out=[]
    xmin=P[:,:,0].min(1); xmax=P[:,:,0].max(1)
    cx=P[:,:,:2].mean(1)
    sel=np.flatnonzero((mat==M.SI)&(cx[:,0]>A.BODY_BOTTOM+1e-12)&(xmin<=x0)&(xmax>=x0))
    for e in sel:
        p=P[e,:3]
        ys=[]
        for k in range(3):
            q1,q2=p[k],p[(k+1)%3]
            if abs(q1[0]-q2[0])<1e-15: continue
            t=(x0-q1[0])/(q2[0]-q1[0])
            if 0<=t<=1: ys.append(q1[1]+t*(q2[1]-q1[1]))
        if len(ys)<2: continue
        y0,y1=min(ys),max(ys)
        if y1-y0<1e-12: continue
        yy=np.linspace(y0,y1,Ns)
        vals=np.array([u_at(p,u[tris[e]],x0,y)[0] for y in yy])
        ind=vals<=eps
        i=0
        while i<Ns:
            if not ind[i]: i+=1; continue
            j=i
            while j+1<Ns and ind[j+1]: j+=1
            lo=y0 if i==0 else yy[i-1]; hi=yy[i]
            for _ in range(60):
                mid=0.5*(lo+hi)
                if u_at(p,u[tris[e]],x0,mid)[0]<=eps: hi=mid
                else: lo=mid
            ys_t=0.5*(lo+hi)
            lo=yy[j]; hi=y1 if j==Ns-1 else yy[j+1]
            for _ in range(60):
                mid=0.5*(lo+hi)
                if u_at(p,u[tris[e]],x0,mid)[0]<=eps: lo=mid
                else: hi=mid
            ye_t=0.5*(lo+hi)
            out.append((ys_t,ye_t)); i=j+1
    out.sort()
    merged=[]
    for a_,b_ in out:
        if merged and a_-merged[-1][1]<2e-3: merged[-1]=(merged[-1][0],max(merged[-1][1],b_))
        else: merged.append((a_,b_))
    return merged

for nm,f in (("L1b","p2_L1b_scaled.npz"),("L1","p2_L1_scaled.npz"),("mid","p2_p3_rmid_polished.npz")):
    z=np.load(f)
    for x0 in (2.35,2.41):
        old=clean(cut_intervals(z["nodes"],z["tris"],z["mat"],z["u"],x0,1e-8,True))[0]
        new=corrected(z["nodes"],z["tris"],z["mat"],z["u"],x0,1e-8)
        print(f"{nm:4s} x={x0}: OLD {[(round(a,4),round(b,4)) for a,b in old]}")
        print(f"          NEW {[(round(a,4),round(b,4)) for a,b in new]}")
