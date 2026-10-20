import json
d = json.load(open("p2_size_audit.json"))
for m, v in d.items():
    a = v.get("A_core"); b = v.get("B_core")
    print(f"{m:4s} A_core med={a['med']:.5f}  B_core med={b['med']:.5f}")
