from pathlib import Path
import pandas as pd

CANDIDATES = Path(
    "results/reports/prioritized_disagreement_review.tsv"
)
SUBMISSIONS = Path(
    "data/processed/hboc_clinvar_submissions.tsv"
)
OUTPUT = Path(
    "results/reports/evidence_provenance_report.tsv"
)

def unique_values(series):
    values = sorted({
        str(value).strip()
        for value in series
        if str(value).strip() not in {"", "-", "NA", "nan"}
    })
    return " | ".join(values)

def main():
    for path in (CANDIDATES, SUBMISSIONS):
        if not path.exists():
            raise FileNotFoundError(
                f"Required file missing: {path}\n"
                "Run the previous project steps first."
            )

    candidates = pd.read_csv(
        CANDIDATES, sep="\t", dtype=str, keep_default_na=False
    )
    submissions = pd.read_csv(
        SUBMISSIONS, sep="\t", dtype=str, keep_default_na=False
    )

    required = {
        "VariationID", "MappedGeneSymbol", "SCV", "Submitter",
        "ClinicalSignificance", "ReviewStatus",
        "DateLastEvaluated", "ReportedPhenotypeInfo"
    }
    missing = required - set(submissions.columns)
    if missing:
        raise ValueError(
            f"Submission file missing columns: {sorted(missing)}"
        )

    candidate_ids = set(candidates["VariationID"])
    relevant = submissions[
        submissions["VariationID"].isin(candidate_ids)
    ].copy()

    # Summarize submission-level provenance for each variant.
    summaries = []
    for variation_id, group in relevant.groupby("VariationID"):
        summaries.append({
            "VariationID": variation_id,
            "SubmissionRows": len(group),
            "DistinctSCVs": group["SCV"].replace(
                {"": pd.NA, "-": pd.NA}
            ).nunique(),
            "DistinctSubmitters": group["Submitter"].replace(
                {"": pd.NA, "-": pd.NA}
            ).nunique(),
            "SubmissionClassifications": unique_values(
                group["ClinicalSignificance"]
            ),
            "SubmissionReviewStatuses": unique_values(
                group["ReviewStatus"]
            ),
            "SubmissionEvaluationDates": unique_values(
                group["DateLastEvaluated"]
            ),
            "SubmissionConditions": unique_values(
                group["ReportedPhenotypeInfo"]
            ),
        })

    provenance = pd.DataFrame(summaries)

    # Avoid collisions with similarly named candidate-report columns.
    provenance = provenance.rename(columns={
        "SubmissionRows": "ProvenanceSubmissionRows",
        "DistinctSCVs": "ProvenanceDistinctSCVs",
        "DistinctSubmitters": "ProvenanceDistinctSubmitters",
    })

    report = candidates.merge(
        provenance, on="VariationID", how="left", validate="many_to_one"
    )

    for column in (
        "ProvenanceSubmissionRows",
        "ProvenanceDistinctSCVs",
        "ProvenanceDistinctSubmitters",
    ):
        report[column] = report[column].fillna(0).astype(int)

    report.to_csv(OUTPUT, sep="\t", index=False)

    print("Evidence provenance report created.")
    print(f"Candidate groups: {len(report):,}")
    print(f"Candidates with matching submissions: "
          f"{report['ProvenanceSubmissionRows'].gt(0).sum():,}")
    print(f"Output: {OUTPUT}")

if __name__ == "__main__":
    main()
