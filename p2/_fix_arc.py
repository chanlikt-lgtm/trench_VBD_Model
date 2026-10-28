p="p2_arc_region.py"; s=open(p,encoding="utf-8").read()
s=s.replace("    if area < AREA_TOL or depth >= MAX_DEPTH:", "    if area < AB.AREA_TOL or depth >= MAX_DEPTH:")
s=s.replace('    out = dict(region="R_arc: |dist to canonical arc| <= 0.01 um",',
            '    AB.AREA_TOL = 1e-13\n    out = dict(region="R_arc: |dist to canonical arc| <= 0.01 um",')
s=s.replace('AREA_TOL=AREA_TOL, levels={})', 'AREA_TOL=AB.AREA_TOL, levels={})')
open(p,"w",encoding="utf-8").write(s)
print("patched to honour AB.AREA_TOL")
