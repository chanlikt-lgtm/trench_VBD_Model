import numpy as np
import p2_front_conditioning as F
z=np.load("p2_p3_rmid_polished.npz")
nodes,tris,u=z["nodes"],z["tris"],z["u"]
print("u stats: min %.3e max %.3e" % (u.min(), u.max()))
print("phi stats: min %.3e" % z["phi"].min())
order,_=F.detect_order(nodes,tris)
pr=F.profile(nodes,tris,u,2.35,order)
ys,us,ds=pr
print("profile y range %.4f..%.4f  n=%d" % (ys.min(),ys.max(),len(ys)))
print("u along cut: min %.3e max %.3e" % (us.min(), us.max()))
for k in range(0,len(ys),max(1,len(ys)//24)):
    print("  y=%+.4f  u=%.6e  du/dy=%+.3e" % (ys[k],us[k],ds[k]))
