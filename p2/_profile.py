import numpy as np
def prof(path, label):
    z = np.load(path); nodes, tris = z["nodes"], z["tris"]
    p = nodes[tris][:, :3]
    det = ((p[:,1,0]-p[:,0,0])*(p[:,2,1]-p[:,0,1]) -
           (p[:,2,0]-p[:,0,0])*(p[:,1,1]-p[:,0,1]))
    h = np.sqrt(4*(0.5*np.abs(det))/np.sqrt(3))
    cx = p[:,:,0].mean(1); cy = p[:,:,1].mean(1)
    print(f"\n=== {label} : median h in x-bins, |y|<=0.30 ===")
    m0 = np.abs(cy) <= 0.30
    edges = np.arange(0.0, 3.01, 0.1)
    row = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = m0 & (cx >= a) & (cx < b)
        row.append(f"{(0.5*(a+b)):.2f}:{np.median(h[m]):.4f}" if m.sum() > 20 else f"{0.5*(a+b):.2f}:--")
    for i in range(0, len(row), 6):
        print("  " + "  ".join(row[i:i+6]))
prof("p2_rmid_p2mesh.npz", "mid")
prof("p2_L1_p2mesh.npz", "L1")
prof("p2_L1b_p2mesh.npz", "L1b")
