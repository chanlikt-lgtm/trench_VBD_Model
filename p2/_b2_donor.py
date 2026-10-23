import numpy as np
import p2_exact_donor as D
import json
z=np.load("p2_B2_p2mesh.npz")
err,tot,f=D.test(z["nodes"],z["tris"],z["mat"],"B2")
m=json.load(open("p2_B2_mesh.json"))
print(f"exact charged (mesh json) = {m['exact_charged']:.12f}")
print(f"sum f / rho               = {tot:.12f}")
print(f"match: {abs(tot-m['exact_charged'])<1e-12}")
