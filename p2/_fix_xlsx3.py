p="p2_make_xlsx.py"; s=open(p,encoding="utf-8").read()
s=s.replace("d01['dA01']", "first")
open(p,"w",encoding="utf-8").write(s)
print("replaced remaining d01['dA01'] references")
