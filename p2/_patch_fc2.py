p="p2_front_conditioning.py"
s=open(p,encoding="utf-8").read()

old_start = s.index("            per_eps = {}")
old_end = s.index("        # endpoints, gradients, conditioning")
new = '''            per_eps = {}
            for eps in EPSS:
                iv2 = cut_intervals(z["nodes"], z["tris"], z["mat"], z["u"], x0, eps, True)
                m2, s2 = clean(iv2)
                per_eps[f"{eps:g}"] = [[float(a), float(b)] for a, b in m2]
            rec[name] = dict(extract_1e8=[[float(a), float(b)] for a, b in main],
                             poly=per_eps, short=int(len(short)))
            print(f"\\n{name}: extractor eps=1e-8 {[(round(a,4),round(b,4)) for a,b in main]}")
            for eps in EPSS[1:]:
                print(f"   eps={eps:g}: {[(round(a,4),round(b,4)) for a,b in per_eps[f'{eps:g}']]}")

'''
s = s[:old_start] + new + s[old_end:]

s = s.replace('''        print(f"{'mesh':5s} {'endpoint':>9s} {'u':>10s} {'du/dy':>10s} {'|dy| for du=1e-8':>17s}")
        for name, _ in MESHES:
            pr = prof[name]
            for (a, b) in rec[name]["poly"]["1e-08"]:
                for y0 in (a, b):
                    g, uu = grad_at(pr[0], pr[1], pr[2], y0, 1e-8)
                    amb = 1e-8/abs(g) if abs(g) > 1e-14 else float("inf")
                    print(f"{name:5s} {y0:+9.4f} {uu:10.2e} {g:+10.3e} {amb:17.3e}")
                    rec[name].setdefault("cond", []).append(
                        dict(y=y0, du=uu, dudy=g, dy_for_1e8=amb))''',
'''        print(f"{'mesh':5s} {'endpoint':>9s} {'u(y*)':>11s} {'du/dy':>11s} {'dy per 1e-8':>13s}")
        for name, _ in MESHES:
            pr = prof[name]
            ys, us, ds = pr
            for (a, b) in rec[name]["poly"]["1e-08"]:
                for y0 in (a, b):
                    uu = float(np.interp(y0, ys, us))
                    g = float(np.interp(y0, ys, ds))
                    amb = 1e-8/abs(g) if abs(g) > 1e-12 else float("inf")
                    print(f"{name:5s} {y0:+9.4f} {uu:11.3e} {g:+11.3e} {amb:13.3e}")
                    rec[name].setdefault("cond", []).append(
                        dict(y=y0, u=uu, dudy=g, dy_for_1e8=amb))''')
open(p,"w",encoding="utf-8").write(s)
print("patched")
