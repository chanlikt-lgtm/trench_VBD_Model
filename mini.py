import gmsh, time
t0=time.time()
gmsh.initialize()
print("initialized", round(time.time()-t0,2), flush=True)
gmsh.model.add("t")
occ=gmsh.model.occ
occ.addRectangle(0,0,0,2.42,1.35)
occ.synchronize()
gmsh.model.mesh.generate(2)
n,_ ,_= gmsh.model.mesh.getNodes()
print("nodes", len(n), "total", round(time.time()-t0,2), flush=True)
gmsh.finalize()
print("done", flush=True)
