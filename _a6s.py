import json
d=json.load(open("p2_3_spatial_S0/CP-L2-A6_arc_region.json"))
print("AREA_TOL", d["AREA_TOL"], "| region", d["region"])
for k,v in d["levels"].items():
    print(f"{k}: arc [{v['arc_lower']:.12e}, {v['arc_upper']:.12e}] w {v['arc_width']:.2e} cand {v['candidates']} "
          f"| comp [{v['comp_lower']:.12f}, {v['comp_upper']:.12f}]")
