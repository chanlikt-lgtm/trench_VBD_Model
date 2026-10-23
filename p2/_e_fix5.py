import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import p2_plot_E as P
nodes,tris,mat,u,Emv = P.load()
x,y=nodes[:,0],nodes[:,1]; tv=tris[:,:3]
fig,ax=plt.subplots(figsize=(5.2,6.4))
keep=(x[tv].mean(1)>=2.20)&(x[tv].mean(1)<=2.43)&(np.abs(y[tv].mean(1))<=0.20)
tp=ax.tripcolor(x,y,tv[keep],facecolors=Emv[keep],shading="flat",cmap="turbo",vmin=0,vmax=0.20)
ax.tricontour(x,y,tv,u,levels=[1e-8],colors="w",linewidths=1.1)
ax.set_xlim(2.20,2.43); ax.set_ylim(-0.20,0.20); ax.set_aspect("equal")
ax.set_xlabel("x (um)"); ax.set_ylabel("y (um)")
ax.set_title("trench interior: |E| 0-0.20 MV/cm\nwhite = neutral boundary")
fig.colorbar(tp,ax=ax,label="|E| (MV/cm)")
fig.tight_layout(); fig.savefig("figE5_interior.png",dpi=170); plt.close(fig); print("figE5_interior.png")
