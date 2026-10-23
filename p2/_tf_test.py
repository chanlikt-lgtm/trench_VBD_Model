import numpy as np, gmsh, sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import mesh_stage2b as M
import p2_audit as A

def build(hg, ht, npts):
    gmsh.initialize(); gmsh.option.setNumber("General.Verbosity", 0)
    arcs=[]
    orig = gmsh.model.occ.addCircleArc
    def patched(*a, **k):
        t = orig(*a, **k); arcs.append(t); return t
    gmsh.model.occ.addCircleArc = patched
    groups = M.build(hg, ht)
    gmsh.model.occ.addCircleArc = orig
    for t in arcs:
        gmsh.model.mesh.setTransfiniteCurve(t, npts)
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.SecondOrderLinear", 1)
    gmsh.option.setNumber("Mesh.MeshSizeMin", ht)
    gmsh.option.setNumber("Mesh.MeshSizeMax", hg)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.model.mesh.generate(2)
    ntags, nc, _ = gmsh.model.mesh.getNodes()
    nodes = np.column_stack([nc[0::3], nc[1::3]])
    t2i = {int(t): i for i,t in enumerate(ntags)}
    acc=[]; mats=[]
    for code, lst in groups.items():
        for face in lst:
            et,_,en = gmsh.model.mesh.getElements(2, face)
            for e,nd in zip(et,en):
                if e==9:
                    acc.append(np.vectorize(t2i.get)(np.asarray(nd,dtype=int).reshape(-1,6)))
                    mats.append(np.full(len(acc[-1]), code))
    gmsh.finalize()
    tris=np.vstack(acc); mat=np.concatenate(mats)
    P=nodes[tris][:,:3]
    det=((P[:,1,0]-P[:,0,0])*(P[:,2,1]-P[:,0,1])-(P[:,2,0]-P[:,0,0])*(P[:,1,1]-P[:,0,1]))
    area=0.5*np.abs(det)
    si=mat==M.SI
    return len(nodes), len(tris), float(area[si].sum()), len(arcs)

for npts in (65,):
    a=build(0.0085,0.0034,npts)
    b=build(0.00601,0.002404,npts)
    print(f"transfinite arcs n={npts}")
    print(f"  L0 nodes={a[0]} elems={a[1]} Si={a[2]:.15f} arcs={a[3]}")
    print(f"  L1 nodes={b[0]} elems={b[1]} Si={b[2]:.15f}")
    print(f"  Si delta = {b[2]-a[2]:+.6e}")
