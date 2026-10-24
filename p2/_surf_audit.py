import gmsh, sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import mesh_stage2b as M
gmsh.initialize(); gmsh.option.setNumber("General.Verbosity", 0)
groups = M.build(0.0085, 0.0034)
surf = [abs(t) for d,t in gmsh.model.getEntities(2)]
used = sorted({abs(f) for lst in groups.values() for f in lst})
print("model 2D surfaces:", len(surf), sorted(surf))
print("extracted faces :", len(used), used)
print("not extracted   :", sorted(set(surf)-set(used)))
for s in sorted(set(surf)-set(used)):
    x0,y0,z0,x1,y1,z1 = gmsh.model.getBoundingBox(2, s)
    print(f"  surface {s}: bbox x[{x0:.3f},{x1:.3f}] y[{y0:.3f},{y1:.3f}]")
gmsh.finalize()
