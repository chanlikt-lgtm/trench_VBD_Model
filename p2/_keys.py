import numpy as np
for f in ("p2_p3_r3_polished.npz","p2_p3_r4_polished.npz","p2_p3_rmid_polished.npz","p2_L1_scaled.npz","p2_L1b_scaled.npz"):
    z=np.load(f)
    print(f, list(z.keys()), "n=", len(z["nodes"]) if "nodes" in z else "?")
