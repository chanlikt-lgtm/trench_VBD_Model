import numpy as np
for f in ("p2_L1_osqp_checkpoint.npz","p2_L1b_osqp_checkpoint.npz"):
    z=np.load(f); y=z["y"]
    print(f, "y min/max/absmax: %+.3e %+.3e %.3e" % (y.min(), y.max(), np.abs(y).max()),
          " nonzero:", int((np.abs(y)>1e-12).sum()))
