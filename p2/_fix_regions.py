p="p2_area_regions.py"; s=open(p,encoding="utf-8").read()
old="""def classify(ph, rect):
    x0, x1, y0, y1 = rect
    xs = ph[:, 0]; ys = ph[:, 1]"""
new="""def classify(ph, rect):
    x0, x1, y0, y1 = rect
    P = np.asarray(ph, dtype=float).reshape(-1, 2)
    xs = P[:, 0]; ys = P[:, 1]"""
assert old in s
s=s.replace(old,new)
open(p,"w",encoding="utf-8").write(s)
print("patched classify")
