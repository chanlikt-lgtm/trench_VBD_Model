import numpy as np
import p2_plot_E as P
nodes,tris,mat,u,Emv = P.load()
x,y=nodes[:,0],nodes[:,1]; tv=tris[:,:3]
cx=x[tv].mean(1); cy=y[tv].mean(1)
idx=np.argsort(Emv)[-8:][::-1]
print("top |E| elements (centroid):")
for k in idx:
    print(f"  {Emv[k]:.3f} MV/cm at x={cx[k]:.4f} y={cy[k]:+.4f}")
si = mat==P.M.SI
top_si = np.argsort(np.where(si,Emv,-1))[-5:][::-1]
print("top within Si:")
for k in top_si:
    print(f"  {Emv[k]:.3f} MV/cm at x={cx[k]:.4f} y={cy[k]:+.4f} mat={mat[k]}")
