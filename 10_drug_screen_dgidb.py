"""
10_drug_screen_dgidb.py - exploratory drug-gene interaction screen (DGIdb v5, GraphQL).

Replaces the old L1000CDS2 / DGIdb v2 REST script (v2 REST is retired).
Input : results/tables/hub_genes.csv (372 ME2 hub genes), results/tables/final_gene_panel.csv
Output: results/tables/dgidb_hub_interactions.csv, results/tables/dgidb_summary.txt
Paper : 372 hubs -> 94 with >=1 interaction, 2,488 rows, 1,950 distinct drugs; none for the 4 panel genes.
This is a lookup of KNOWN interactions - NOT a prediction that any drug treats aneurysm rupture.

NOTE: written against the public DGIdb GraphQL schema (endpoint below) but not executed
in the authoring sandbox (network blocked). If the API returns a schema error, open
https://dgidb.org/api and adjust the field names in QUERY.
Run: pip install requests pandas ; python scripts/10_drug_screen_dgidb.py
"""
import json, time
import pandas as pd
import requests

URL = "https://dgidb.org/api/graphql"
QUERY = """
query GetInteractions($genes: [String!]) {
  genes(names: $genes) {
    nodes {
      name
      interactions {
        drug { name conceptId }
        interactionScore
        interactionTypes { type directionality }
        sources { fullName }
      }
    }
  }
}
"""

def query(genes):
    r = requests.post(URL, json={"query": QUERY, "variables": {"genes": genes}}, timeout=120)
    r.raise_for_status()
    js = r.json()
    if "errors" in js:
        raise RuntimeError(json.dumps(js["errors"])[:500])
    return js["data"]["genes"]["nodes"]

def rows_from(nodes):
    out = []
    for n in nodes:
        for i in n.get("interactions") or []:
            out.append({
                "gene": n["name"],
                "drug": i["drug"]["name"],
                "drug_concept_id": i["drug"].get("conceptId"),
                "interaction_score": i.get("interactionScore"),
                "interaction_types": ";".join(t["type"] for t in (i.get("interactionTypes") or []) if t.get("type")),
                "sources": ";".join(s["fullName"] for s in (i.get("sources") or [])),
            })
    return out

hubs = pd.read_csv("results/tables/hub_genes.csv")["gene"].tolist()
panel = pd.read_csv("results/tables/final_gene_panel.csv")["gene"].tolist()

rows, B = [], 40
for k in range(0, len(hubs), B):
    rows += rows_from(query(hubs[k:k + B]))
    time.sleep(0.5)
df = pd.DataFrame(rows)
df.to_csv("results/tables/dgidb_hub_interactions.csv", index=False)

panel_rows = pd.DataFrame(rows_from(query(panel)))
summary = (
    f"Hub genes queried: {len(hubs)}\n"
    f"Hub genes with >=1 interaction: {df['gene'].nunique() if len(df) else 0}\n"
    f"Interaction rows: {len(df)}\n"
    f"Distinct drugs: {df['drug'].nunique() if len(df) else 0}\n"
    f"Panel genes ({', '.join(panel)}) with >=1 interaction: "
    f"{panel_rows['gene'].nunique() if len(panel_rows) else 0}\n"
)
open("results/tables/dgidb_summary.txt", "w").write(summary)
print(summary)
