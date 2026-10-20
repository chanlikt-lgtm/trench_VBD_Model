p="p2_make_figs.py"
s=open(p,encoding="utf-8").read()
s=s.replace('"increments grow as the mesh is refined; they do not decay"',
            '"increment reduction reverses at the finest local refinement"')
open(p,"w",encoding="utf-8").write(s)
print("caption updated")
