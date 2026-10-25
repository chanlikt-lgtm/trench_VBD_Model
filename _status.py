import json
a=json.load(open("p2_3_spatial_S0/A_S0-L0.json"))
c=json.load(open("p2_3_spatial_S0/CP-M2-C_L0.json"))
print("A_S0-L0:", a["lower"], a["upper"], "width", a["width"])
print("L0 F7:", c["f7"]["stationarity"], c["f7"]["complementarity"], "gate", c["f7"]["gate"],
      "| minres_info", c.get("minres_info"))
print("OSQP:", c["osqp"]["iters"], "iters, peak RSS", c["peak_rss_mb"], "MB")
