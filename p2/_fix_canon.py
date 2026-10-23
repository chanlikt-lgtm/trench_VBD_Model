p="p2_mesh_S0.py"; s=open(p,encoding="utf-8").read()
s=s.replace('DIST_MIN, DIST_MAX = 0.04, 0.06',
            'DIST_MIN, DIST_MAX = 0.04, 0.06\nN_ARC = 65          # canonical chord count per circular arc (64 chords)')
s=s.replace('''    occ = gmsh.model.occ
    tags = [occ.addPoint(float(x), float(y), 0.0, hb) for x, y in sel]
    spl = occ.addSpline(tags)
    occ.synchronize()''',
'''    occ = gmsh.model.occ
    tags = [occ.addPoint(float(x), float(y), 0.0, hb) for x, y in sel]
    spl = occ.addSpline(tags)
    occ.synchronize()

    # Canonical geometry: force a fixed chord count on every circular arc so the
    # polygonal domain is identical at all levels (CP-M1a). Query by curve type
    # because OCC fragmentation renumbers entities.
    arcs = [abs(t) for d, t in gmsh.model.getEntities(1)
            if gmsh.model.getType(1, abs(t)) == "Circle"]
    for t in arcs:
        gmsh.model.mesh.setTransfiniteCurve(t, N_ARC)
    print(f"[{prefix}] canonical arcs: {len(arcs)} x {N_ARC} nodes", flush=True)''')
s=s.replace('rec = dict(level=level, rk=rk,',
            'rec = dict(level=level, rk=rk, n_arc=N_ARC, n_arcs=len(arcs),')
open(p,"w",encoding="utf-8").write(s)
print("patched canonical arcs")
