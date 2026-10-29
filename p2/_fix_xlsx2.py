p="p2_make_xlsx.py"; s=open(p,encoding="utf-8").read()
s=s.replace('d01 = json.load(open(D / "DA_S0_01.json"))\nfirst = d12["dA01"] if "dA01" in d12 else d01["delta"]',
            'd01 = json.load(open(D / "DA_S0_01.json"))')
s=s.replace('d12 = json.load(open(D / "DA_S0_12.json"))',
            'd12 = json.load(open(D / "DA_S0_12.json"))\nfirst = d12["dA01"] if "dA01" in d12 else d01["delta"]')
open(p,"w",encoding="utf-8").write(s)
print("reordered definition")
