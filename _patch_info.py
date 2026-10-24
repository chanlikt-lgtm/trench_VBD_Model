p="p2/p2_cpm2c_l0.py"; s=open(p,encoding="utf-8").read()
s=s.replace("        zt, info = spla.minres(Hop, D*b, x0=z0, rtol=1e-8, maxiter=8000)",
            "        zt, info = spla.minres(Hop, D*b, x0=z0, rtol=1e-8, maxiter=8000)\n        minres_info.append(int(info))")
s=s.replace("    repairs = []\n    best = None",
            "    repairs = []\n    minres_info = []\n    best = None")
s=s.replace('               active=int(len(W)), repairs=repairs,',
            '               active=int(len(W)), repairs=repairs, minres_info=minres_info,')
open(p,"w",encoding="utf-8").write(s)
print("patched future runs to record minres info")
