import numpy as np, json
import p2_area_bounds as AB
for tag,f in (("mid","p2_p3_rmid_polished.npz"),):
    z=np.load(f)
    lo,up,nch,nm=AB.bounds(z["nodes"],z["tris"],z["mat"],z["u"])
    print(f"{tag}: whole-domain lower={lo:.12f} upper={up:.12f} gap={up-lo:.2e} charged={nch} mixed={nm}")
    print(f"  recorded p2_area_bounds.json value: ", json.load(open("p2_area_bounds.json"))["mid"])
