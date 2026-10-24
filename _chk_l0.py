import json
d=json.load(open("p2_3_spatial_S0/CP-M2-C_L0.json"))
print("minres_info:", d.get("minres_info"))
a=json.load(open("p2_3_spatial_S0/A_S0-L0.json"))
print("A_S0-L0:", a["lower"], a["upper"], "width", a["width"], "mixed", a["mixed"], a["seconds"], "s")
