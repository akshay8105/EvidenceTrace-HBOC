#!/usr/bin/env python3

"""Filter ClinVar submissions and identify potential condition-specific disagreements."""

from pathlib import Path
import gzip
import pandas as pd

AGGREGATE = Path("data/processed/hboc_clinvar_grch38.tsv")
SUBMISSIONS = Path("data/raw/submission_summary.txt.gz")

OUT_DATA = Path("data/processed/hboc_clinvar_submissions.tsv")
OUT_REPORTS = Path("results/reports")

TARGET_GENES = {"BRCA1", "BRCA2", "PALB2", "ATM", "CHEK2"}
CHUNK_SIZE = 100_000

COLUMNS = [
    "VariationID",
    "ClinicalSignificance",
    "DateLastEvaluated",
    "Description",
    "SubmittedPhenotypeInfo",
    "ReportedPhenotypeInfo",
    "ReviewStatus",
    "CollectionMethod",
    "OriginCounts",
    "Submitter",
    "SCV",
    "SubmittedGeneSymbol",
    "ExplanationOfInterpretation",
    "ContributesToAggregateClassification",
]


def main():
    if not AGGREGATE.exists() or not SUBMISSIONS.exists():
        raise FileNotFoundError(
            "Required input missing. Check the aggregate and submission files."
        )

    OUT_DATA.parent.mkdir(parents=True, exist_ok=True)
    OUT_REPORTS.mkdir(parents=True, exist_ok=True)
    OUT_DATA.unlink(missing_ok=True)

    aggregate = pd.read_csv(
        AGGREGATE,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        usecols=["VariationID", "GeneSymbol"],
    )

    aggregate = aggregate[
        aggregate["GeneSymbol"].isin(TARGET_GENES)
    ].copy()

    aggregate["VariationID"] = aggregate["VariationID"].str.strip()
    aggregate = aggregate[aggregate["VariationID"] != ""]

    # Preserve all gene mappings in case an identifier maps to
    # more than one target gene in the aggregate file.
    gene_map = (
        aggregate.groupby("VariationID")["GeneSymbol"]
        .agg(lambda x: ";".join(sorted(set(x))))
        .to_dict()
    )
    target_ids = set(gene_map)

    if not target_ids:
        raise RuntimeError("No target VariationIDs found in aggregate dataset.")

    written = False
    total_rows = 0
    selected_rows = 0
    chunks_seen = 0

    # The actual tabular header is line 19 (1-based).
    # Skip the 18 preceding documentation lines.
    with gzip.open(SUBMISSIONS, "rt", encoding="utf-8") as handle:
        reader = pd.read_csv(
            handle,
            sep="\t",
            skiprows=18,
            dtype=str,
            keep_default_na=False,
            chunksize=CHUNK_SIZE,
            low_memory=False,
        )

        for chunk in reader:
            chunks_seen += 1
            chunk.columns = [
                str(c).strip().lstrip("#").strip()
                for c in chunk.columns
            ]
            total_rows += len(chunk)

            missing = set(COLUMNS) - set(chunk.columns)
            if missing:
                raise ValueError(
                    f"Missing expected columns: {sorted(missing)}. "
                    f"Found: {chunk.columns.tolist()}"
                )

            chunk["VariationID"] = chunk["VariationID"].str.strip()
            selected = chunk[
                chunk["VariationID"].isin(target_ids)
            ].copy()

            if not selected.empty:
                selected["MappedGeneSymbol"] = selected["VariationID"].map(
                    gene_map
                )
                selected = selected[
                    COLUMNS + ["MappedGeneSymbol"]
                ]
                selected.to_csv(
                    OUT_DATA,
                    sep="\t",
                    index=False,
                    mode="a",
                    header=not written,
                )
                written = True
                selected_rows += len(selected)

            if chunks_seen % 10 == 0:
                print(
                    f"Scanned {total_rows:,} submissions; "
                    f"matched {selected_rows:,} so far."
                )

    if not written:
        raise RuntimeError(
            "No submission records matched target VariationIDs. "
            "Check the file header and VariationID values."
        )

    print("\nSubmission extraction complete.")
    print(f"Target VariationIDs: {len(target_ids):,}")
    print(f"Submission rows scanned: {total_rows:,}")
    print(f"Matched submission rows: {selected_rows:,}")
    print(f"Output: {OUT_DATA}")

    # Summarize submissions by mapped gene.
    submissions = pd.read_csv(
        OUT_DATA, sep="\t", dtype=str, keep_default_na=False
    )

    # Normally each ID maps to one gene. Keep ambiguous mappings visible.
    gene_counts = []
    for gene in sorted(TARGET_GENES):
        subset = submissions[
            submissions["MappedGeneSymbol"].apply(
                lambda value: gene in value.split(";")
            )
        ]
        gene_counts.append({
            "GeneSymbol": gene,
            "SubmissionRows": len(subset),
            "DistinctVariationIDs": subset["VariationID"].nunique(),
            "DistinctSCVs": subset["SCV"].nunique(),
        })

    pd.DataFrame(gene_counts).to_csv(
        OUT_REPORTS / "submission_counts_by_gene.tsv",
        sep="\t",
        index=False,
    )

    # Classification frequencies by gene and submitter-level label.
    classification_counts = []
    for _, row in pd.DataFrame(gene_counts).iterrows():
        gene = row["GeneSymbol"]
        subset = submissions[
            submissions["MappedGeneSymbol"].apply(
                lambda value: gene in value.split(";")
            )
        ]
        counts = (
            subset.groupby("ClinicalSignificance")
            .size()
            .reset_index(name="SubmissionRows")
        )
        counts.insert(0, "GeneSymbol", gene)
        classification_counts.append(counts)

    pd.concat(classification_counts, ignore_index=True).to_csv(
        OUT_REPORTS / "submission_classifications_by_gene.tsv",
        sep="\t",
        index=False,
    )

    # Look for multiple distinct classifications for the same VariationID
    # and non-missing reported condition among submissions marked "yes".
    contributing = submissions[
        submissions["ContributesToAggregateClassification"]
        .str.strip().str.lower().eq("yes")
    ].copy()

    condition = contributing["ReportedPhenotypeInfo"].str.strip()
    condition_missing = condition.str.lower().isin(
        {"", "na", "-", "not provided", "not specified"}
    )
    contributing = contributing[~condition_missing].copy()

    group_columns = ["VariationID", "ReportedPhenotypeInfo"]
    disagreement_rows = []

    for (variation_id, phenotype), group in contributing.groupby(
        group_columns, dropna=False
    ):
        labels = sorted({
            label.strip()
            for label in group["ClinicalSignificance"]
            if label.strip() and label.strip().lower() not in {"na", "-"}
        })

        if len(labels) > 1:
            disagreement_rows.append({
                "VariationID": variation_id,
                "MappedGeneSymbol": ";".join(
                    sorted(set(gene_map.get(variation_id, "").split(";")))
                ),
                "ReportedPhenotypeInfo": phenotype,
                "DistinctClassificationCount": len(labels),
                "Classifications": " | ".join(labels),
                "ContributingSubmissionRows": len(group),
                "DistinctSCVs": group["SCV"].nunique(),
                "Submitters": " | ".join(sorted({
                    x.strip() for x in group["Submitter"]
                    if x.strip() and x.strip() != "-"
                })),
            })

    disagreements = pd.DataFrame(
        disagreement_rows,
        columns=[
            "VariationID", "MappedGeneSymbol", "ReportedPhenotypeInfo",
            "DistinctClassificationCount", "Classifications",
            "ContributingSubmissionRows", "DistinctSCVs", "Submitters",
        ],
    )

    disagreements.to_csv(
        OUT_REPORTS / "potential_condition_specific_disagreements.tsv",
        sep="\t",
        index=False,
    )

    print("\nSubmission counts by gene:")
    print(pd.DataFrame(gene_counts).to_string(index=False))
    print(
        "\nPotential disagreement groups with multiple classifications: "
        f"{len(disagreements):,}"
    )
    print(
        "Note: these are screening candidates, not independently "
        "validated clinical conflicts."
    )
    print("\nReports written to:")
    print(f"  {OUT_REPORTS / 'submission_counts_by_gene.tsv'}")
    print(f"  {OUT_REPORTS / 'submission_classifications_by_gene.tsv'}")
    print(
        f"  {OUT_REPORTS / 'potential_condition_specific_disagreements.tsv'}"
    )


if __name__ == "__main__":
    main()
