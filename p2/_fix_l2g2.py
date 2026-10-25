p="p2_cpl2_geom.py"; s=open(p,encoding="utf-8").read()

old = '''    bn2 = np.unique(np.sort(np.vstack([
        t2[m2 == M.SI][:, [0, 1]], t2[m2 == M.SI][:, [1, 2]],
        t2[m2 == M.SI][:, [2, 0]]]), axis=1)[:, :].ravel())
    pts = n2[np.unique(bn2)]'''
new = '''    s2 = t2[m2 == M.SI][:, :3]
    e2 = np.sort(np.vstack([s2[:, [0, 1]], s2[:, [1, 2]], s2[:, [2, 0]]]), axis=1)
    u2, c2 = np.unique(e2, axis=0, return_counts=True)
    pts = n2[np.unique(u2[c2 == 1].ravel())]'''
assert old in s
s = s.replace(old, new)

# arc-adjacent region: elements whose centroid lies close to the frozen arc circles
old2 = '''    defs = {"band<=0.04": d <= 0.04,'''
new2 = '''    # distance to the two canonical arc circles (radius R about the arc centres)
    R = abs(M.A_BOT - M.A_ARC_END)
    darc = np.minimum(np.abs(np.hypot(cx - M.TIP_START, cy - M.A_ARC_END) - R),
                      np.abs(np.hypot(cx - M.TIP_START, cy + M.A_ARC_END) - R))
    defs = {"arc<=0.01": darc <= 0.01,
            "band<=0.04": d <= 0.04,'''
assert old2 in s
s = s.replace(old2, new2)
open(p,"w",encoding="utf-8").write(s)
print("patched boundary filter and arc region")
