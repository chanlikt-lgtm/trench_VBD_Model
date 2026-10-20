import numpy as np
import p2_audit as A

M=A.M
# --- synthetic check: u = 2x - 3y on a unit triangle -------------------
P=np.array([[0.,0.],[1.,0.],[0.,1.]])
def uv(x,y): return 2.0*x-3.0*y
nodes6=np.array([[0,0],[1,0],[0,1],[.5,0],[.5,.5],[0,.5]],float)
ue=np.array([uv(*p) for p in nodes6])
def ev(P,ue,x0,y):
    Mx=np.array([[P[1,0]-P[0,0],P[2,0]-P[0,0]],[P[1,1]-P[0,1],P[2,1]-P[0,1]]])
    Mi=np.linalg.inv(Mx)
    st=Mi@(np.array([x0,y])-P[0]); l1,l2=st; l3=1-l1-l2
    return float(A.shape_vals(l1,l2)@ue)
h=1e-6
fd=(ev(P,ue,0.25,0.25+h)-ev(P,ue,0.25,0.25-h))/(2*h)
print(f"synthetic u=2x-3y: exact du/dy=-3, finite-diff={fd:.9f}")

# --- real meshes: containing element + finite-difference gradient ------
def eval_pt(nodes,tris,mat,u,x0,y):
    P=nodes[tris]
    xmin=P[:,:,0].min(1); xmax=P[:,:,0].max(1)
    sel=np.flatnonzero((xmin<=x0+1e-12)&(xmax>=x0-1e-12)&(mat==M.SI))
    outs=[]
    for e in sel:
        p=P[e,:3]
        Mx=np.array([[p[1,0]-p[0,0],p[2,0]-p[0,0]],[p[1,1]-p[0,1],p[2,1]-p[0,1]]])
        try: Mi=np.linalg.inv(Mx)
        except np.linalg.LinAlgError: continue
        st=Mi@(np.array([x0,y])-p[0]); l1,l2=st; l3=1-l1-l2
        if l1>=-1e-12 and l2>=-1e-12 and l3>=-1e-12:
            ue=u[tris[e]]
            v=lambda yy: float(A.shape_vals(*np.linalg.solve(Mx,np.array([x0,yy])-p[0]))@ue)
            outs.append((v(y),(v(y+h)-v(y-h))/(2*h)))
    return outs[0] if outs else None

for name,f in (("L1b","p2_L1b_scaled.npz"),("L1","p2_L1_scaled.npz")):
    z=np.load(f); nodes,tris,mat,u=z["nodes"],z["tris"],z["mat"],z["u"]
    print(f"\n=== {name} x=2.35 ===")
    for y in (-0.100,-0.090,-0.085,-0.080,-0.0774,-0.075,-0.0718,-0.068,-0.060):
        r=eval_pt(nodes,tris,mat,u,2.35,y)
        print(f"  y={y:+.4f}  u={r[0]:+.6e}  du/dy={r[1]:+.4e}" if r else f"  y={y:+.4f} outside Si")
