import numpy as np
import p2_audit as A
M=A.M; xb=A.BODY_BOTTOM
for tag,f in (("B0","p2_B0_p2mesh.npz"),("B1","p2_B1_p2mesh.npz")):
    z=np.load(f); nodes,tris,mat=z["nodes"],z["tris"],z["mat"]
    P=nodes[tris][:,:3]; cx=P[:,:,0].mean(1); si=mat==M.SI
    xmax=P[:,:,0].max(1); xmin=P[:,:,0].min(1)
    charged=si&(cx>xb)                 # centroid rule: constraint set + old load
    contrib=si&(xmax>xb)               # exact-load candidate set
    gained=contrib&~charged            # centroid said no, exact says yes
    lost=charged&~contrib              # impossible
    print(f"{tag}: charged(constraints)={charged.sum()} contrib(exact)={contrib.sum()} "
          f"gained={gained.sum()} lost={lost.sum()}")
