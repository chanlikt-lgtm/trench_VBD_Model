import numpy as np
import p2_L1b_polish as P

d = P.build()
chk = np.load("p2_L1b_osqp_checkpoint.npz")
u = chk["x"]; y = chk["y"]
C_F, cD = d["C_F"], d["cD"]
slack = C_F @ u + cD
print(f"rows={len(slack)}  min slack={slack.min():+.3e}  max={slack.max():+.3e}")
for thr in (1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7):
    print(f"  slack<{thr:.0e}: {(slack<thr).sum():7d}   & lam>1e-10: {((slack<thr)&(y>1e-10)).sum():7d}")
print("lam counts:")
for thr in (1e-12, 1e-10, 1e-8, 1e-6, 1e-4):
    print(f"  lam>{thr:.0e}: {(y>thr).sum():7d}")
W1 = np.flatnonzero(slack < 1e-7)
W2 = np.flatnonzero((slack < 1e-9) & (y > 1e-10))
W3 = np.flatnonzero((slack < 1e-8) & (y > 1e-10))
print(f"\n|W slack<1e-7|={len(W1)}  |intersection 1e-9|={len(W2)}  |intersection 1e-8|={len(W3)}")
np.save("_W_slack.npy", W1); np.save("_W_int.npy", W2)
