p="p2_cpm2_assembly.py"; s=open(p,encoding="utf-8").read()
s=s.replace('''    fin = dict(K=np.isfinite(K.data).all(), f=np.isfinite(f).all(),
               C=np.isfinite(C.data).all(), Kff=np.isfinite(Kff.data).all())''',
'''    fin = dict(K=bool(np.isfinite(K.data).all()), f=bool(np.isfinite(f).all()),
               C=bool(np.isfinite(C.data).all()), Kff=bool(np.isfinite(Kff.data).all()))''')
open(p,"w",encoding="utf-8").write(s)
print("patched bool serialisation")
