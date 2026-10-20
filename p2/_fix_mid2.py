p="p2_area_regions.py"; s=open(p,encoding="utf-8").read()
s=s.replace("""def mid3(P, Q):
    return ((P[0]+Q[0])/2, (P[1]+Q[1])/2, (P[2]+Q[2])/2)""",
"""def mid3(P, Q):
    return ((P[0]+Q[0])/2, (P[1]+Q[1])/2, (P[2]+Q[2])/2)


def mid2(P, Q):
    return ((P[0]+Q[0])/2, (P[1]+Q[1])/2)""")
s=s.replace("m01 = mid3(ph[0], ph[1]); m12 = mid3(ph[1], ph[2]); m20 = mid3(ph[2], ph[0])",
            "m01 = mid2(ph[0], ph[1]); m12 = mid2(ph[1], ph[2]); m20 = mid2(ph[2], ph[0])")
open(p,"w",encoding="utf-8").write(s)
print("patched mid2")
