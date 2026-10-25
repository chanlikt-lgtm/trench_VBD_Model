import json
d=json.load(open("p2_3_spatial_S0/CP-L1_L1.json"))
print("status", d["status"], "| peak RSS", d["peak_rss_mb"], "MB")
print("minres_info", d["minres_info"], "| repairs", d["repairs"])
print("min avail", min(v["avail"] for v in d["marks"].values()), "MB")
print("F7", d["f7"])
