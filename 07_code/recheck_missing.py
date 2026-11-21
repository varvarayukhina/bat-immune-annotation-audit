import json, sys, time, urllib.request, urllib.error
TAX = {"Mmusculus": 10090, "Raegyptiacus": 9407}
BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
n_found = 0
for line in open(sys.argv[1]).read().splitlines()[1:]:
    t, s, why = line.split("\t")
    url = f"{BASE}/gene/symbol/{s}/taxon/{TAX[t]}?page_size=1000"
    time.sleep(0.25)
    try:
        d = json.loads(urllib.request.urlopen(url, timeout=60).read())
    except urllib.error.HTTPError as e:
        print(t, s, "HTTP", e.code); continue
    reps = d.get("reports", [])
    if reps:
        n_found += 1
    print(t, s, d.get("total_count"), [(r["gene"]["symbol"], r["gene"]["gene_id"]) for r in reps])
print("rechecked symbols that now return something:", n_found)
