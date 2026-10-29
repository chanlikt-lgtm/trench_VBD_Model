p="p2_make_xlsx.py"; s=open(p,encoding="utf-8").read()
s=s.replace('d01 = json.load(open(D / "DA_S0_01.json"))',
            'd01 = json.load(open(D / "DA_S0_01.json"))\nfirst = d12["dA01"] if "dA01" in d12 else d01["delta"]')
s=s.replace('d01["dA01"]','first')
open(p,"w",encoding="utf-8").write(s)
print("patched increment source")
