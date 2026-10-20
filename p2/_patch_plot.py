p="p2_plot_trench.py"
s=open(p,encoding="utf-8").read()
s=s.replace('    return nodes[used], o2n[t], phi[used], u[used]',
            '    return nodes[used], o2n[t][:, :3], phi[used], u[used]')
s=s.replace('        n, t, phi, u = load(f)\n        tri = mtri.Triangulation(n[:, 0], n[:, 1], t)',
            '        n, t, phi, u = load(f)\n        tri = mtri.Triangulation(n[:, 0], n[:, 1], t)')
open(p,"w",encoding="utf-8").write(s)
print("patched")
