"""
09_drug_repurposing.py  (Role 5 owns this)
Queries L1000CDS2 with the up/down hub-gene signature to find compounds
predicted to REVERSE the rupture-associated expression pattern, then
cross-checks the top candidates against DGIdb for known druggability.

Needs Stage 4's hub_genes.csv / module_trait_relationships.csv - the
deconvolution half of Stage 8 can run without this, but this script can't.
"""

import requests
import pandas as pd
import os

os.makedirs("results/tables", exist_ok=True)

hub_df = pd.read_csv("results/tables/hub_genes.csv")
module_trait = pd.read_csv("results/tables/module_trait_relationships.csv")
top_module_corr_sign = module_trait.iloc[0]["correlation"]

if top_module_corr_sign > 0:
    up_genes = hub_df[hub_df["module_membership"] > 0]["gene"].tolist()
    down_genes = hub_df[hub_df["module_membership"] < 0]["gene"].tolist()
else:
    down_genes = hub_df[hub_df["module_membership"] > 0]["gene"].tolist()
    up_genes = hub_df[hub_df["module_membership"] < 0]["gene"].tolist()

print(f"Up-in-rupture genes: {len(up_genes)}")
print(f"Down-in-rupture genes: {len(down_genes)}")

if len(up_genes) < 5 or len(down_genes) < 5:
    print("WARNING: L1000CDS2 performs best with ~10-150 genes per direction.")
    print("Consider using module_gene_full_table.csv instead if either list is too short.")

L1000CDS2_URL = "https://maayanlab.cloud/L1000CDS2/query"
payload = {
    "data": {"upGenes": up_genes, "dnGenes": down_genes},
    "config": {"aggravate": False, "searchMethod": "geneSet", "share": False,
               "combination": False, "db-version": "latest"}
}

print("\nQuerying L1000CDS2 (reversal mode)...")
drug_df = pd.DataFrame()
try:
    resp = requests.post(L1000CDS2_URL, json=payload, timeout=60)
    resp.raise_for_status()
    result = resp.json()
    drug_df = pd.DataFrame(result.get("topMeta", []))
    if not drug_df.empty:
        drug_df.to_csv("results/tables/l1000cds2_candidates.csv", index=False)
        print(f"Saved {len(drug_df)} candidate compounds.")
        print(drug_df.head(10))
    else:
        print("L1000CDS2 returned no results - check gene symbols are valid HGNC.")
except requests.exceptions.RequestException as e:
    print(f"L1000CDS2 request failed: {e}")

if not drug_df.empty and "pert_desc" in drug_df.columns:
    top_drugs = drug_df["pert_desc"].dropna().unique().tolist()[:20]
    DGIDB_URL = "https://dgidb.org/api/v2/interactions.json"
    dgidb_rows = []
    print(f"\nCross-checking top {len(top_drugs)} candidates against DGIdb...")
    for drug in top_drugs:
        try:
            r = requests.get(DGIDB_URL, params={"drugs": drug}, timeout=20)
            r.raise_for_status()
            for m in r.json().get("matchedTerms", []):
                for interaction in m.get("interactions", []):
                    dgidb_rows.append({
                        "drug": drug, "gene": interaction.get("geneName"),
                        "interaction_types": ",".join(interaction.get("interactionTypes", [])),
                        "sources": ",".join(interaction.get("sources", []))
                    })
        except requests.exceptions.RequestException as e:
            print(f"  DGIdb lookup failed for {drug}: {e}")
    dgidb_df = pd.DataFrame(dgidb_rows)
    dgidb_df.to_csv("results/tables/dgidb_crosscheck.csv", index=False)
    print(f"Saved DGIdb cross-check ({len(dgidb_df)} rows).")
    print("NOTE: DGIdb's API endpoint has changed across versions - if this")
    print("404s, check https://www.dgidb.org/api for the current endpoint.")
else:
    print("\nSkipping DGIdb cross-check - no L1000CDS2 candidates to check.")

print("\n=== Stage 9 complete ===")
print("Molecular docking against the top hub protein is an optional stretch")
print("goal, not included here - do this manually if time allows.")
