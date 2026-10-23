import numpy as np
import p2_plot_E as P
nodes,tris,mat,u,Emv = P.load()
x,y = nodes[:,0], nodes[:,1]; tv=tris[:,:3]
cx=x[tv].mean(1); cy=y[tv].mean(1)
m=(cx>=2.20)&(cx<=2.50)&(np.abs(cy)<=0.20)
e=Emv[m]
print(f"window elements {m.sum()}")
for q in (50,75,90,95,98,99,100):
    print(f"  p{q}: {np.percentile(e,q):.3f} MV/cm")
