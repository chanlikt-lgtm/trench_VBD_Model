p="p2_cpl2_geom.py"; s=open(p,encoding="utf-8").read()
s=s.replace('    for k in data["L2"]["regions"]:\n        if k in data["L1"]["regions"]:\n            a = data["L1"]["regions"][k]; b = data["L2"]["regions"][k]',
            '    for k in out["L2"]["regions"]:\n        if k in out["L1"]["regions"]:\n            a = out["L1"]["regions"][k]; b = out["L2"]["regions"][k]')
open(p,"w",encoding="utf-8").write(s)
print("patched regions reference")
