import numpy as np
import p2_audit as A

# unit right triangle: v0=(0,0), v1=(1,0), v2=(0,1)
P = np.array([[0.,0.],[1.,0.],[0.,1.]])
V = np.array([[0.,0.],[1.,0.],[0.,1.],[.5,0.],[.5,.5],[0.,.5]])
u = 2.0*V[:,0] - 3.0*V[:,1]          # exact nodal values of 2x-3y
Mx = np.array([[P[1,0]-P[0,0], P[2,0]-P[0,0]],
               [P[1,1]-P[0,1], P[2,1]-P[0,1]]])
Mi = np.linalg.inv(Mx)

def conv_A(x,y):
    a,b = Mi @ (np.array([x,y]) - P[0])
    return float(A.shape_vals(a,b) @ u)          # current code path

def conv_B(x,y):
    a,b = Mi @ (np.array([x,y]) - P[0])
    L1 = 1.0-a-b; L2 = a; L3 = b                # barycentric of v0,v1,v2
    N = np.array([L1*(2*L1-1), L2*(2*L2-1), L3*(2*L3-1),
                  4*L1*L2, 4*L2*L3, 4*L3*L1])
    return float(N @ u)

print("point        exact        conv_A       conv_B")
for (x,y) in ((0.25,0.25),(0.6,0.2),(0.1,0.7),(0.45,0.45)):
    print(f"({x:.2f},{y:.2f})   {2*x-3*y:+.6f}   {conv_A(x,y):+.6f}   {conv_B(x,y):+.6f}")

h=1e-6
for name,fn in (("A",conv_A),("B",conv_B)):
    fd=(fn(0.25,0.25+h)-fn(0.25,0.25-h))/(2*h)
    print(f"conv_{name}: du/dy finite-diff = {fd:+.6f}   (exact -3)")
