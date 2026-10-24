import json, hashlib, platform, sys
import numpy, scipy
try:
    import osqp; OSQP_V = getattr(osqp,"__version__","?")
except Exception as e:
    OSQP_V = f"missing:{e}"
try:
    import gmsh; GMSH_V = getattr(gmsh,"__version__","?")
except Exception as e:
    GMSH_V = f"missing:{e}"

def h(p):
    return hashlib.sha256(open(p,"rb").read()).hexdigest()

mesh = {t: h(f"../p2_3_spatial_S0/levels/S0_{t}_p2mesh.npz") for t in ("L0","L1")}
code = {f: h(f) for f in ("p2_audit.py","p2_exact_donor.py","p2_area_bounds.py",
                          "p2_exact_metrics.py","p2_nested_consistency.py",
                          "p2_bernstein.py","p2_mesh_S0.py")}
band = h("p2_band_contour.npz")

man = dict(
  checkpoint="CP-M2-A",
  status="FROZEN",
  provenance=dict(design="e937926", repair="1c8d4df", preflight="82a4a1b/fefbba1"),
  mesh=dict(L0=dict(sha256=mesh["L0"], nodes=232552, compacted=231274,
                    elements=115123, charged=59322, bernstein_rows=355932),
            L1=dict(sha256=mesh["L1"], nodes=460914, compacted=459288,
                    elements=228917, charged=118290, bernstein_rows=709740)),
  band_sha256=band,
  geometry=dict(policy="G-fixed polygonal domain",
                arcs=4, nodes_per_arc=65, chords_per_90deg=64,
                boundary="frozen L1b contour, 625 resampled points"),
  discretization=dict(space="conforming straight-sided P2 (6-node)",
                      element_order=2, second_order_linear=1),
  donor=dict(rule="exact clipped-element P2 integration, T cap {x>=0.54}",
             quadrature="3-point degree-2 in parent barycentrics",
             centroid_fallback=False),
  obstacle=dict(domain="charged silicon only", rows="6 per charged element",
                construction="Bernstein nonnegativity, unchanged"),
  solver=dict(osqp=dict(eps_abs=1e-8, eps_rel=1e-8, max_iter=30000,
                        scaling=10, polishing=False, verbose=False),
              polish=dict(method="matrix-free scaled MINRES + residual correction",
                          rtol=1e-8, maxiter=8000, nref=6,
                          correction_sign="[-stat, +slackW]")),
  f7=dict(stationarity=1e-8, complementarity=1e-10, min_slack=-1e-9,
          violations=0, lambda_min=-1e-9, oxide=1e-9, overshoot=1e-9),
  checkpointing=dict(save=("x","y","idx"), reload_verify=True,
                     before="polishing", after="OSQP status=solved"),
  software=dict(python=platform.python_version(), numpy=numpy.__version__,
                scipy=scipy.__version__, osqp=OSQP_V, gmsh=GMSH_V,
                platform=platform.platform()),
  code_sha256=code,
  deviations_from_validated_X=["mesh: canonical 64-chord arcs (geometry freeze)",
                               "DOF count differs from B-X family",
                               "no other change to formulation or tolerances"],
)
json.dump(man, open("../p2_3_spatial_S0/CP-M2-A_manifest.json","w"), indent=2)
print("software:", man["software"])
print("mesh sha L0:", mesh["L0"][:32])
print("mesh sha L1:", mesh["L1"][:32])
print("wrote CP-M2-A_manifest.json")
