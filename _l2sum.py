import json
d=json.load(open("p2_3_spatial_S0/CP-L2_L2.json"))
m=d["marks"]
print("runtime %.1f min | peak %.0f MB | min avail %d MB" %
      (m["after_polish"]["t"]/60, d["peak_rss_mb"], min(v["avail"] for v in m.values())))
print("OSQP %d iters | active %d | repairs %s | minres %s" %
      (d["osqp"]["iters"], d["active"], d["repairs"], d["minres_info"]))
