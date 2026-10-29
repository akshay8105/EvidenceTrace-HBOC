#!/usr/bin/env python3

"""Generate an exploratory evidence audit of filtered ClinVar records."""

from pathlib import Path
import pandas as pd

INPUT = Path("data/processed/hboc_clinvar_grch38.tsv")
OUT = Path("results/reports")
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT, sep="\t", dtype=str, keep_default_na=False)

required = {
    "GeneSymbol", "VariationID", "ClinicalSignificance",
    "ReviewStatus", "NumberSubmitters", "Assembly",
}
missing = required.difference(df.columns)
if missing:
    raise ValueError(f"Missing required columns: {sorted(missing)}")

# Detect explicit aggregate conflict labels; do not infer conflicts
# merely because different conditions appear in the same row.
df["HasClassificationConflict"] = (
    df["ClinicalSignificance"]
    .str.contains("conflicting", case=False, na=False)
)

# Summary 1: rows and distinct VariationIDs by gene.
gene_summary = (
    df.groupby("GeneSymbol", dropna=False)
    .agg(
        RecordRows=("VariationID", "size"),
        DistinctVariationIDs=("VariationID", "nunique"),
        ExplicitConflictRows=("HasClassificationConflict", "sum"),
    )
    .reset_index()
)
gene_summary.to_csv(OUT / "clinvar_gene_audit.tsv", sep="\t", index=False)

# Summary 2: distribution of aggregate clinical-significance labels.
significance_summary = (
    df.groupby(["GeneSymbol", "ClinicalSignificance"], dropna=False)
    .size()
    .reset_index(name="RecordRows")
    .sort_values(["GeneSymbol", "RecordRows"], ascending=[True, False])
)
significance_summary.to_csv(
    OUT / "clinvar_significance_audit.tsv", sep="\t", index=False
)

# Summary 3: review-status distribution.
review_summary = (
    df.groupby(["GeneSymbol", "ReviewStatus"], dropna=False)
    .size()
    .reset_index(name="RecordRows")
    .sort_values(["GeneSymbol", "RecordRows"], ascending=[True, False])
)
review_summary.to_csv(
    OUT / "clinvar_review_status_audit.tsv", sep="\t", index=False
)

# Summary 4: examples of records explicitly labelled as conflicting.
conflicts = df.loc[df["HasClassificationConflict"]].copy()
conflicts = conflicts[
    [
        "GeneSymbol", "VariationID", "Name", "ClinicalSignificance",
        "ReviewStatus", "NumberSubmitters", "LastEvaluated",
        "RCVaccession", "PhenotypeList",
    ]
]
conflicts.to_csv(
    OUT / "clinvar_conflict_examples.tsv", sep="\t", index=False
)

print("ClinVar audit complete.")
print(f"Total rows: {len(df):,}")
print(f"Unique VariationIDs across file: {df['VariationID'].nunique():,}")
print(f"Explicit conflict-labelled rows: {int(df['HasClassificationConflict'].sum()):,}")
print("\nGene summary:")
print(gene_summary.to_string(index=False))
print("\nTop 15 significance labels:")
print(
    significance_summary.groupby("ClinicalSignificance")["RecordRows"]
    .sum().sort_values(ascending=False).head(15).to_string()
)
print("\nReports written to:")
for path in sorted(OUT.glob("clinvar_*_audit.tsv")):
    print(f"  {path}")
print(f"  {OUT / 'clinvar_conflict_examples.tsv'}")
