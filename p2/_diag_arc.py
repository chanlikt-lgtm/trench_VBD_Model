import gmsh, sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import mesh_stage2b as M

arcs=[]
orig = gmsh.model.occ.addCircleArc
def patched(*a, **k):
    t = orig(*a, **k); arcs.append(t); return t
gmsh.model.occ.addCircleArc = patched
gmsh.initialize(); gmsh.option.setNumber("General.Verbosity", 0)
groups = M.build(0.0085, 0.0034)
gmsh.model.occ.addCircleArc = orig
curves = [abs(t) for d,t in gmsh.model.getEntities(1)]
print("captured arc tags:", arcs)
print("model curve count:", len(curves), "first few:", sorted(curves)[:8])
print("all captured present:", all(t in curves for t in arcs))
print("model name:", gmsh.model.getCurrent())
gmsh.finalize()
