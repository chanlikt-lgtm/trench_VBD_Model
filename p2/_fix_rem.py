p="p2_area_regions.py"; s=open(p,encoding="utf-8").read()
old = """        def rec(pc, bar, ph, ar, depth):
            if any(classify(ph, r) > 0 for r in REGIONS.values()):
                return
            sc = sub_coeff(pc, bar[0], bar[1], bar[2])
            if all(ci <= 0.0 for ci in sc):
                acc[0] += ar; return
            if all(ci > 0.0 for ci in sc):
                return
            if ar < AREA_TOL or depth >= MAX_DEPTH:
                acc[1] += ar; return"""
new = """        def rec(pc, bar, ph, ar, depth):
            ca = classify(ph, REGIONS["A"]); cb = classify(ph, REGIONS["B"])
            if ca > 0 or cb > 0:
                return                      # belongs to A or B
            straddle = (ca == 0) or (cb == 0)
            sc = sub_coeff(pc, bar[0], bar[1], bar[2])
            if not straddle:
                if all(ci <= 0.0 for ci in sc):
                    acc[0] += ar; return
                if all(ci > 0.0 for ci in sc):
                    return
            if ar < AREA_TOL or depth >= MAX_DEPTH:
                acc[1] += ar; return"""
assert old in s
s = s.replace(old, new)
open(p,"w",encoding="utf-8").write(s)
print("patched remainder straddle handling")
