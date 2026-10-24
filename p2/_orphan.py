import numpy as np
from scipy.spatial import cKDTree
band=np.load("p2_band_contour.npz")["pts"]
tree=cKDTree(band)
for t in ("L0","L1"):
    z=np.load(f"../p2_3_spatial_S0/levels/S0_{t}_p2mesh.npz")
    nodes,tris=z["nodes"],z["tris"]
    used=np.unique(tris.ravel()); orb=np.setdiff1d(np.arange(len(nodes)),used)
    d,_=tree.query(nodes[orb])
    print(f"{t}: orphans {len(orb)}  dist to band: max {d.max():.2e} p50 {np.median(d):.2e}")
    print(f"     x range {nodes[orb,0].min():.3f}..{nodes[orb,0].max():.3f}  "
          f"y range {nodes[orb,1].min():.3f}..{nodes[orb,1].max():.3f}")
