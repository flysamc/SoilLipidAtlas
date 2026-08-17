#!/usr/bin/env python3
"""Supplementary Figure 7 panel b - per-phylum evidence tiers recomputed on the
LOCKED release ncbi-phylum-2026-08-04-v1.

Scope decision (author-directed, 2026-08-11): panel b describes SAMPLING DEPTH,
i.e. what was profiled, so it is NOT restricted to the 16 analysis units. It
shows the full collection carrying the corrected NCBI phylum assignments used
to run the strict analysis - the 19 collection phyla.

What the corrected assignment changes relative to the published panel
(all verified to reconcile sample-for-sample against the legacy table):
  Amoebozoa 16 + Heterolobosea 1  -> Discosea 7 + Evosea 5 + Heterolobosea 5
  Euryarchaeota 14 + Methanobacteriota 2 -> Methanobacteriota 16
  Crenarchaeota 1 + Thermoproteota 2     -> Thermoproteota 3
  Bryophyta 4 + Marchantiophyta 3 + Trachaeophyta 5 + Magnoliophyta 1
    + Charophyta 1                       -> Streptophyta 14
  Mucoromycota 8                         -> Mucoromycota 7 + Mortierellomycota 1

Counting basis: taxonomy_scope == 'core_candidate' (samples with a valid NCBI
phylum ancestor), 168 POS and 195 NEG. Descriptive-only labels (Bicosoecida,
Rootnodules, Mixed) have no phylum-rank ancestor and are reported as a footnote
row rather than as phyla; viral units are out of study scope (author
instruction, 2026-08-03). The published panel counted 169 POS samples because it
displayed Bicosoecida as a unit; 168 + 1 = 169.

Tier rule, as recovered from the analysis provenance notes:
  A robust      >= 5 samples, >= 3 genera, >= 2 batches
  B moderate    >= 3 samples
  C preliminary  = 2 samples
  D anecdotal    = 1 sample
Applied to positive-mode depth, consistently. The published tier column does NOT
follow this rule for three units (Euryarchaeota, Pseudomonadota, Mucoromycota
were published as B while meeting every Tier A threshold) and no further
criterion separating them could be recovered; every such deviation is written to
tier_changes_vs_published.csv so it is visible rather than silent.

Arthropoda keeps its published asterisk: it meets the Tier A sample and species
thresholds but derives from a single acquisition batch.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RELEASE_ROOT = PROJECT_ROOT / "outputs" / "analysis" / "ncbi-phylum-2026-08-04-v1"
TAX = RELEASE_ROOT / "taxonomy"
OUTPUT = RELEASE_ROOT / "suppfig7_evidence_tiers_release_2026-08-11_v1"
PUBLISHED = (PROJECT_ROOT / "manuscript_2_clean" / "06_figures" / "figures_r"
             / "data" / "supp_annotation" / "evidence_tiers.csv")
POLICY = PROJECT_ROOT / "paper2_repro" / "config" / "taxonomy_policy.json"

SINGLE_BATCH_ASTERISK = "Arthropoda"


def tier_of(samples: int, genera: int, batches: int) -> str:
    if samples >= 5 and genera >= 3 and batches >= 2:
        return "A"
    if samples >= 3:
        return "B"
    if samples == 2:
        return "C"
    return "D"


def main() -> int:
    if OUTPUT.exists():
        sys.exit(f"Refusing to overwrite {OUTPUT} - delete or rename it first.")

    summary = json.loads((TAX / "taxonomy_summary.json").read_text())
    policy = json.loads(POLICY.read_text(encoding="utf-8"))
    assert policy["status"] == "locked"
    collection = sorted(summary["collection_phyla"])
    assert len(collection) == 19

    pos = pd.read_csv(TAX / "sample_metadata_POS_ncbi_phylum.csv")
    neg = pd.read_csv(TAX / "sample_metadata_NEG_ncbi_phylum.csv")
    pos_core = pos[pos["taxonomy_scope"] == "core_candidate"]
    neg_core = neg[neg["taxonomy_scope"] == "core_candidate"]
    assert len(pos_core) == 168 and len(neg_core) == 195

    # Taxonomic breadth column.  The published `pos_species` column is NOT
    # reproducible and is internally inconsistent: Euryarchaeota was published
    # as 4 species across 9 distinct genera (impossible), and Mucoromycota's 8
    # equals its raw sample count rather than its 6 distinct species (the
    # duplicate arises from a case typo, "Rhizopus Oligosporus" vs
    # "oligosporus").  The `species` field is also null for 42 of 168 core
    # samples and carries the literal value "unclassified", while `genus`
    # sometimes holds a bare genus and sometimes a binomial plus strain.
    # Distinct GENUS is the finest rank the metadata determines for every
    # sample, so breadth is reported as genera and the tier rule uses it.
    pos_core = pos_core.copy()
    pos_core["genus_norm"] = (pos_core["genus"].astype(str).str.strip()
                              .str.lower().str.split().str[0])
    p = pos_core.groupby("ncbi_phylum").agg(
        pos_samples=("original_header", "nunique"),
        pos_genera=("genus_norm", "nunique"),
        pos_batches=("batch", "nunique"))
    n = neg_core.groupby("ncbi_phylum").agg(
        neg_samples=("sample_col", "nunique"),
        neg_batches=("batch", "nunique"))

    ev = (p.join(n, how="outer").reindex(collection).fillna(0).astype(int)
          .reset_index().rename(columns={"ncbi_phylum": "phylum"}))
    ev["tier"] = [tier_of(r.pos_samples, r.pos_genera, r.pos_batches)
                  for r in ev.itertuples()]
    ev["note"] = ["*" if r.phylum == SINGLE_BATCH_ASTERISK and r.pos_batches == 1
                  else "" for r in ev.itertuples()]
    ev["ecological_group"] = ev["phylum"].map(policy["ecological_group"])
    ev["in_strict_analysis"] = ev["phylum"].isin(summary["analysis_phyla"]).astype(int)
    ev = ev.sort_values(["tier", "pos_samples", "phylum"],
                        ascending=[True, False, True]).reset_index(drop=True)

    assert ev["pos_samples"].sum() == 168, ev["pos_samples"].sum()
    assert ev["neg_samples"].sum() == 195, ev["neg_samples"].sum()

    pub = pd.read_csv(PUBLISHED)
    merged = ev.merge(pub[["phylum", "tier"]], on="phylum", how="left",
                      suffixes=("", "_published"))
    changes = merged[merged["tier"] != merged["tier_published"]][
        ["phylum", "pos_samples", "pos_genera", "pos_batches",
         "tier", "tier_published"]].copy()
    changes["reason"] = [
        "unit did not exist under the published labels"
        if pd.isna(r.tier_published) else
        "published tier does not follow the documented rule; rule applied consistently"
        for r in changes.itertuples()]

    descriptive = (pos[pos["taxonomy_scope"] == "descriptive_only"]
                   .groupby("source_phylum").size().rename("pos_samples")
                   .reset_index().rename(columns={"source_phylum": "label"}))
    descriptive["reason"] = descriptive["label"].map(
        policy["descriptive_only_labels"]).fillna("descriptive-only label")

    OUTPUT.mkdir(parents=True)
    ev.to_csv(OUTPUT / "evidence_tiers_release.csv", index=False)
    changes.to_csv(OUTPUT / "tier_changes_vs_published.csv", index=False)
    descriptive.to_csv(OUTPUT / "descriptive_only_labels.csv", index=False)
    (OUTPUT / "RUN_SUMMARY.json").write_text(json.dumps({
        "status": "recomputed from locked release metadata (pure counting); "
                  "published tier column is not reproducible for 3 units",
        "taxonomy_release": summary["taxonomy_release"],
        "scope": "collection phyla (sampling depth), NOT restricted to the 16 "
                 "analysis units - author direction 2026-08-11",
        "n_collection_phyla": len(collection),
        "n_analysis_phyla": len(summary["analysis_phyla"]),
        "counting_basis": "taxonomy_scope == core_candidate",
        "n_pos_samples": int(ev["pos_samples"].sum()),
        "n_neg_samples": int(ev["neg_samples"].sum()),
        "published_pos_sample_total": 169,
        "sample_total_reconciliation":
            "168 core-candidate POS samples + 1 Bicosoecida sample shown as a "
            "unit in the published panel = 169",
        "tier_counts": ev["tier"].value_counts().sort_index().to_dict(),
        "published_tier_counts": pub["tier"].value_counts().sort_index().to_dict(),
        "n_tier_changes_vs_published": int(len(changes)),
        "tier_rule": {"A": ">=5 samples, >=3 genera, >=2 batches",
                      "B": ">=3 samples", "C": "2 samples", "D": "1 sample"},
        "unreproducible_published_tiers": [
            "Euryarchaeota", "Pseudomonadota", "Mucoromycota"],
        "descriptive_only_excluded_from_phylum_rows": descriptive["label"].tolist(),
    }, indent=2) + "\n", encoding="utf-8")

    print(ev.to_string(index=False))
    print(f"\ntier counts: {ev['tier'].value_counts().sort_index().to_dict()} "
          f"(published {pub['tier'].value_counts().sort_index().to_dict()})")
    print(f"\n{len(changes)} tier changes vs published:")
    print(changes.to_string(index=False))
    print(f"\nWrote {OUTPUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
