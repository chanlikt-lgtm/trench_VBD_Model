p="p2_cpm2_pre.py"; s=open(p,encoding="utf-8").read()

# exact point-to-segment distance instead of sampled points
s=s.replace('''    seg_pts = []
    for a, b in bnd:
        pa, pb = n0[a], n0[b]
        for s in np.linspace(0, 1, 5):
            seg_pts.append(pa + s*(pb-pa))
    tree = cKDTree(np.array(seg_pts))
    d, _ = tree.query(n1[bn1])''',
'''    Aseg = n0[bnd[:, 0]]; Bseg = n0[bnd[:, 1]]
    pts = n1[bn1]
    best = np.full(len(pts), np.inf)
    for i in range(len(Aseg)):
        a = Aseg[i]; ab = Bseg[i]-a; L2 = float(ab@ab)
        t = np.clip(((pts-a) @ ab)/L2, 0.0, 1.0)
        proj = a + t[:, None]*ab
        d = np.linalg.norm(pts-proj, axis=1)
        best = np.minimum(best, d)
    d = best''')

# orientation per material, not a global inversion count
s=s.replace('''    return dict(min_angle=float(ang.min()),''',
'''    code = None
    return dict(min_angle=float(ang.min()),''')
s=s.replace('''def quality(nodes, tris):''','''def orientation(nodes, tris, mat):
    P = nodes[tris][:, :3]
    det = ((P[:, 1, 0]-P[:, 0, 0])*(P[:, 2, 1]-P[:, 0, 1]) -
           (P[:, 2, 0]-P[:, 0, 0])*(P[:, 1, 1]-P[:, 0, 1]))
    out = {}
    for code in np.unique(mat):
        m = mat == code
        out[str(int(code))] = dict(n=int(m.sum()),
                                   pos=int((det[m] > 0).sum()),
                                   neg=int((det[m] < 0).sum()),
                                   zero=int((det[m] == 0).sum()))
    return out


def quality(nodes, tris):''')
s=s.replace('''            quality=quality(nodes, tris),''','''            quality=quality(nodes, tris),
            p2_dofs_compacted=int(len(used)),
            orientation=orientation(nodes, tris, mat),''')
s=s.replace('''            p2_dofs=int(len(nodes)),''','''            p2_dofs_raw=int(len(nodes)),''')
open(p,"w",encoding="utf-8").write(s)
print("patched: exact segment distance, per-material orientation, compacted DOFs")
