import p2_mesh_S0 as S
rec = S.build(2, "../p2_3_spatial_S0/levels")
print("level", rec["level"], "rk", rec["rk"], "band", rec["band"], "bg", rec["background"], "tip", rec["tip"])
