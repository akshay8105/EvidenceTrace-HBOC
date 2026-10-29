from pathlib import Path
import pandas as pd

REPORT = Path("results/reports/evidence_provenance_report.tsv")

def main():
    if not REPORT.exists():
        raise FileNotFoundError(f"Missing report: {REPORT}")

    df = pd.read_csv(REPORT, sep="\t", dtype=str, keep_default_na=False)

    required = {
        "VariationID",
        "MappedGeneSymbol",
        "ReportedPhenotypeInfo",
        "NormalizedClassifications",
        "ReviewPriority",
        "ProvenanceSubmissionRows",
        "ProvenanceDistinctSCVs",
        "ProvenanceDistinctSubmitters",
    }

    missing = required - set(df.columns)
    assert not missing, f"Missing columns: {sorted(missing)}"
    assert len(df) > 0, "Report is empty"
    assert df["VariationID"].str.strip().ne("").all(), (
        "Some rows have missing VariationIDs"
    )

    valid_priorities = {
        "High review priority",
        "Moderate review priority",
        "Routine follow-up",
    }
    unexpected = set(df["ReviewPriority"]) - valid_priorities
    assert not unexpected, f"Unexpected priorities: {unexpected}"

    for column in (
        "ProvenanceSubmissionRows",
        "ProvenanceDistinctSCVs",
        "ProvenanceDistinctSubmitters",
    ):
        values = pd.to_numeric(df[column], errors="raise")
        assert values.ge(0).all(), f"Negative counts in {column}"

    duplicate_keys = df.duplicated(
        subset=["VariationID", "ReportedPhenotypeInfo"]
    ).sum()

    print("Provenance report validation passed.")
    print(f"Rows checked: {len(df):,}")
    print(f"Columns checked: {len(df.columns)}")
    print(f"Duplicate variant-condition rows: {duplicate_keys:,}")
    print("\nReview-priority distribution:")
    print(df["ReviewPriority"].value_counts().to_string())

    if duplicate_keys:
        print(
            "\nNote: duplicate variant-condition rows exist; "
            "inspect them before downstream analyses."
        )

if __name__ == "__main__":
    main()
