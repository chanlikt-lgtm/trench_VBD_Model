import json
p="p2_3_spatial_S0/CP-M2-C_L0.json"
d=json.load(open(p))
d["minres_info"]="NOT_RECORDED"
d["minres_info_note"]=("The MINRES return flag was not serialised in this run. "
                       "F7 numeric certification stands independently; "
                       "the solver-return-code provenance is incomplete. "
                       "Serialisation is mandatory for future solves.")
json.dump(d, open(p,"w"), indent=2)
print("minres_info set to NOT_RECORDED")
