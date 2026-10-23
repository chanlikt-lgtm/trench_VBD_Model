import json, os, hashlib, subprocess
base="p2_3_spatial_S0"
os.makedirs(f"{base}/checkpoints", exist_ok=True)
def h(p):
    return hashlib.sha256(open(p,"rb").read()).hexdigest()
def hh(p): return h(p)[:16].upper()

code={f: hh(f"p2/{f}") for f in ("p2_audit.py","p2_exact_donor.py","p2_area_bounds.py",
      "p2_exact_metrics.py","p2_nested_consistency.py","p2_topology_clean.py","p2_mesh_B2.py")}
code["p2_band_contour.npz"]=hh("p2/p2_band_contour.npz")

man=dict(
  study="P2.3-S0 controlled spatial convergence",
  phase="design_only",
  git_commit=subprocess.check_output(["git","rev-parse","HEAD"]).decode().strip(),
  code_hashes=code,
  protocol_sha256=h(f"{base}/protocol/P2_3_SPATIAL_CONVERGENCE_PROTOCOL_S0.md"),
  model=dict(pde="2D electrostatics, conforming P2", V_D=34.274208536,
             N_D_cm3=6.9e16, rho=1248.5638, body_cut_um=0.54,
             oxide_eps="A.EPS_OX", si_eps="A.EPS_SI"),
  size_field=dict(rule="h_k = r^k h0, r = 2^-1/2, global",
                  band_core_um=0.0040, band_dist_um=0.04, transition_um=0.02,
                  background_um=0.0085, tip_um=0.0034,
                  band_component=dict(n=1, length_um=1.2492,
                                      x=[2.333,2.420], y=[-0.575,0.575],
                                      points=625, hash=code["p2_band_contour.npz"])),
  levels=[dict(level="S0-L0",factor=1.0,band=0.00400,background=0.00850),
          dict(level="S0-L1",factor=0.7071,band=0.0028284,background=0.00601),
          dict(level="S0-L2",factor=0.5,band=0.00200,background=0.00425)],
  solver=dict(osqp="checkpoint-first, eps 1e-8/1e-8, max_iter 30000",
              polish="matrix-free scaled MINRES + residual correction, sign [-stat,+slack]",
              f7="full-set gate, tolerances unchanged"),
  area=dict(estimator="Bernstein branch-and-bound, outward intervals",
            AREA_TOL=1e-13, MAX_DEPTH=24, EPS=1e-8),
  authorizations=dict(mesh_generation=False, solve=False, production=False),
)
json.dump(man, open(f"{base}/protocol/manifests.json","w"), indent=2)

def ck(name,status,phase,evidence,next_action):
    return dict(checkpoint=name,status=status,phase=phase,utc_timestamp=None,
                git_commit=man["git_commit"],source_hashes=code,config_hashes={},
                commands=[],tests_and_metrics=evidence,output_artifacts=[],
                deviations=[],reviewer=None,next_authorized_action=next_action)

tests=json.load(open("p2/p2_s0_tests.json"))
ckpts={
 "CP-S0.0": ck("CP-S0.0","PASS","design_only",
   dict(commits_verified=["15ec811","2478f82","4e5e6d4"],
        chain="CLOSED", solve_running=False), "review CP-S0.1 decisions"),
 "CP-S0.1": ck("CP-S0.1","PASS","design_only",
   dict(target="Bernstein-bounded area of {u<=1e-8} in charged Si",
        geometry="G-fixed polygonal domain",
        geometry_budget="declared, pending arc study",
        front="separate diagnostic series"), "proceed to CP-S0.2"),
 "CP-S0.2": ck("CP-S0.2","PASS","design_only",
   dict(tests=tests, fixture_defects_found=2),
   "review master field (CP-S0.3)"),
 "CP-S0.3": ck("CP-S0.3","PASS_WITH_ESTIMATES","design_only",
   dict(levels=man["levels"],
        resource_anchor="B2-X 125579 elems / 252290 DOFs / 680 MB / 5450 iters / 697 s area",
        estimates=dict(L0="115k elems, ~10 min",L1="230k, ~25 min",L2="460k, ~60 min"),
        pilot_required=True), "assemble S0 package (CP-S0.4)"),
}
for k,v in ckpts.items():
    json.dump(v, open(f"{base}/checkpoints/{k}.json","w"), indent=2)
print("wrote manifests.json and 4 checkpoint records")
print("levels:", [l["level"] for l in man["levels"]])
