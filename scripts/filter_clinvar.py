#!/usr/bin/env python3

"""Filter ClinVar variant summary records for selected HBOC-related genes."""

import gzip
from pathlib import Path

import pandas as pd

INPUT = Path("data/raw/variant_summary.txt.gz")
OUTPUT = Path("data/processed/hboc_clinvar_grch38.tsv")
SUMMARY = Path("results/reports/clinvar_filter_summary.tsv")

TARGET_GENES = {"BRCA1", "BRCA2", "PALB2", "ATM", "CHEK2"}
ASSEMBLY = "GRCh38"
CHUNK_SIZE = 100_000

KEEP_COLUMNS = [
    "#AlleleID",
    "Type",
    "Name",
    "GeneID",
    "GeneSymbol",
    "ClinicalSignificance",
    "ClinSigSimple",
    "LastEvaluated",
    "RS# (dbSNP)",
    "RCVaccession",
    "PhenotypeList",
    "Origin",
    "OriginSimple",
    "Assembly",
    "Chromosome",
    "Start",
    "Stop",
    "ReferenceAllele",
    "AlternateAllele",
    "ReviewStatus",
    "NumberSubmitters",
    "Guidelines",
    "VariationID",
    "PositionVCF",
    "ReferenceAlleleVCF",
    "AlternateAlleleVCF",
    "SCVsForAggregateGermlineClassification",
]


def main():
    if not INPUT.exists():
        raise FileNotFoundError(
            f"Input file not found: {INPUT}. "
            "Download variant_summary.txt.gz first."
        )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)

    # Remove an old output so reruns don't append duplicate data.
    OUTPUT.unlink(missing_ok=True)

    gene_counts = {gene: 0 for gene in sorted(TARGET_GENES)}
    total_rows = 0
    selected_rows = 0
    output_written = False

    with gzip.open(INPUT, "rt", encoding="utf-8") as handle:
        reader = pd.read_csv(
            handle,
            sep="\t",
            dtype=str,
            chunksize=CHUNK_SIZE,
            low_memory=False,
        )

        for chunk_number, chunk in enumerate(reader, start=1):
            total_rows += len(chunk)

            # Some ClinVar files use a leading # in the first column name.
            if "GeneSymbol" not in chunk.columns:
                raise ValueError(
                    "GeneSymbol column not found. "
                    f"Observed columns: {chunk.columns.tolist()[:10]}"
                )

            required = set(KEEP_COLUMNS)
            missing = required.difference(chunk.columns)
            if missing:
                raise ValueError(
                    f"Required columns missing: {sorted(missing)}"
                )

            mask = (
                chunk["GeneSymbol"].isin(TARGET_GENES)
                & chunk["Assembly"].eq(ASSEMBLY)
            )
            selected = chunk.loc[mask, KEEP_COLUMNS].copy()

            if not selected.empty:
                selected_rows += len(selected)

                for gene, count in selected["GeneSymbol"].value_counts().items():
                    gene_counts[gene] += int(count)

                selected.to_csv(
                    OUTPUT,
                    sep="\t",
                    index=False,
                    mode="a",
                    header=not output_written,
                )
                output_written = True

            if chunk_number % 10 == 0:
                print(
                    f"Processed {total_rows:,} records; "
                    f"selected {selected_rows:,} so far."
                )

    if not output_written:
        raise RuntimeError(
            "No records matched the target genes and GRCh38. "
            "Check the input file and its column names."
        )

    # Keep the raw selected records intact. Report possible duplicate
    # representations separately instead of silently deleting them.
    summary = pd.DataFrame(
        [{"GeneSymbol": gene, "SelectedRows": count}
         for gene, count in gene_counts.items()]
    )
    summary.loc[len(summary)] = {
        "GeneSymbol": "TOTAL",
        "SelectedRows": selected_rows,
    }
    summary.to_csv(SUMMARY, sep="\t", index=False)

    print("\nFiltering complete.")
    print(f"Input records processed: {total_rows:,}")
    print(f"Selected GRCh38 rows: {selected_rows:,}")
    print(f"Output: {OUTPUT}")
    print(f"Summary: {SUMMARY}")
    print("\nCounts by gene:")
    for gene, count in gene_counts.items():
        print(f"  {gene}: {count:,}")


if __name__ == "__main__":
    main()
