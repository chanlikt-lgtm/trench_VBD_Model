p="p2_mesh_S0.py"; s=open(p,encoding="utf-8").read()
s=s.replace('out = "p2_3_spatial_S0/levels"','out = os.path.join("..","p2_3_spatial_S0","levels")')
s=s.replace('band = np.load("p2/p2_band_contour.npz")["pts"] if os.path.exists("p2/p2_band_contour.npz") \\\n        else np.load("p2_band_contour.npz")["pts"]',
            'band = np.load("p2_band_contour.npz")["pts"]')
open(p,"w",encoding="utf-8").write(s)
print("patched paths")
