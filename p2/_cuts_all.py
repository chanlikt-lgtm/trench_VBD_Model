import numpy as np
from p2_exact_metrics import cut_intervals
from p2_topology_clean import clean
M=[("r3","p2_p3_r3_polished.npz"),("r4","p2_p3_r4_polished.npz"),
   ("mid","p2_p3_rmid_polished.npz"),("L1b","p2_L1b_scaled.npz"),("L1","p2_L1_scaled.npz")]
for x0 in (2.35, 2.41):
    print(f"--- x={x0} eps=1e-8 ---")
    for n,f in M:
        z=np.load(f)
        iv=cut_intervals(z["nodes"],z["tris"],z["mat"],z["u"],x0,1e-8,True)
        main,short=clean(iv)
        print(f"  {n:4s} " + " ".join(f"({a:+.4f},{b:+.4f})" for a,b in main) + f"  short={len(short)}")
