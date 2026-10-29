from pathlib import Path
import pandas as pd

INPUT = Path("data/processed/hboc_clinvar_submissions.tsv")
OUT = Path("results/reports")
OUT.mkdir(parents=True, exist_ok=True)

GENERIC_CONDITIONS = {
    "", "na", "-", "not provided", "not specified",
    "unknown", "unspecified",
}

def normalize_classification(value):
    """Normalize common labels without replacing the original label."""
    label = str(value).strip().lower()
    label = " ".join(label.split())

    if not label or label in {"na", "-", "not provided"}:
        return "Missing/unspecified"

    if "conflicting classifications" in label:
        return "Conflicting aggregate label"
    if label == "pathogenic/likely pathogenic":
        return "Pathogenic/Likely pathogenic"
    if label == "benign/likely benign":
        return "Benign/Likely benign"
    if label == "likely pathogenic":
        return "Likely pathogenic"
    if label == "pathogenic":
        return "Pathogenic"
    if label == "likely benign":
        return "Likely benign"
    if label == "benign":
        return "Benign"
    if "uncertain significance" in label:
        return "Uncertain significance"
    if "risk factor" in label:
        return "Risk factor"
    if "drug response" in label:
        return "Drug response"
    if "protective" in label:
        return "Protective"
    if "association" in label:
        return "Association"
    if label == "other":
        return "Other"

    return "Other/unmapped: " + label


def main():
    if not INPUT.exists():
        raise FileNotFoundError(f"Missing input: {INPUT}")

    df = pd.read_csv(INPUT, sep="\t", dtype=str, keep_default_na=False)

    required = {
        "VariationID", "ClinicalSignificance", "ReportedPhenotypeInfo",
        "ReviewStatus", "Submitter", "SCV",
        "ContributesToAggregateClassification", "DateLastEvaluated",
        "MappedGeneSymbol",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")

    df["NormalizedClassification"] = df["ClinicalSignificance"].apply(
        normalize_classification
    )

    df["Contributes"] = (
        df["ContributesToAggregateClassification"]
        .str.strip().str.lower().eq("yes")
    )

    condition = df["ReportedPhenotypeInfo"].str.strip()
    condition_key = condition.str.lower()
    df["HasSpecificCondition"] = ~condition_key.isin(GENERIC_CONDITIONS)

    eligible = df[df["Contributes"] & df["HasSpecificCondition"]].copy()

    group_cols = ["VariationID", "ReportedPhenotypeInfo"]
    rows = []

    for (vid, phenotype), group in eligible.groupby(group_cols, dropna=False):
        labels = sorted(set(group["NormalizedClassification"]))
        # Keep only groups with at least two distinct normalized labels.
        if len(labels) < 2:
            continue

        rows.append({
            "VariationID": vid,
            "MappedGeneSymbol": ";".join(sorted(set(
                ";".join(group["MappedGeneSymbol"]).split(";")
            ) - {""})),
            "ReportedPhenotypeInfo": phenotype,
            "NormalizedClassifications": " | ".join(labels),
            "OriginalClassifications": " | ".join(sorted(
                set(group["ClinicalSignificance"])
            )),
            "DistinctClassificationCount": len(labels),
            "ContributingSubmissionRows": len(group),
            "DistinctSCVs": group["SCV"].nunique(),
            "DistinctSubmitters": group["Submitter"].replace(
                {"": pd.NA, "-": pd.NA}
            ).nunique(),
            "Submitters": " | ".join(sorted(set(
                x for x in group["Submitter"] if x and x != "-"
            ))),
            "SCVs": " | ".join(sorted(set(
                x for x in group["SCV"] if x and x != "-"
            ))),
            "ReviewStatuses": " | ".join(sorted(set(
                x for x in group["ReviewStatus"] if x and x != "-"
            ))),
            "EvaluationDates": " | ".join(sorted(set(
                x for x in group["DateLastEvaluated"] if x and x != "-"
            ))),
            "ReasonFlagged": (
                "Multiple normalized classifications among submissions "
                "marked as contributing, for the same VariationID and "
                "reported condition"
            ),
        })

    result = pd.DataFrame(rows)

    output = OUT / "refined_disagreement_candidates.tsv"
    result.to_csv(output, sep="\t", index=False)

    print("Refined disagreement analysis complete.")
    print(f"Eligible contributing submissions: {len(eligible):,}")
    print(f"Candidate groups: {len(result):,}")
    print(f"Report: {output}")

    if not result.empty:
        print("\nCandidate groups by gene:")
        print(result["MappedGeneSymbol"].value_counts().to_string())
        print("\nMost common normalized classifications:")
        print(
            result["NormalizedClassifications"]
            .value_counts().head(10).to_string()
        )
    else:
        print("No groups met the refined screening rule.")

if __name__ == "__main__":
    main()
