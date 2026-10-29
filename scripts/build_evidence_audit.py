from pathlib import Path
import pandas as pd

INPUT = Path("data/processed/hboc_clinvar_grch38.tsv")
OUTPUT = Path("results/reports/evidencetrace_evidence_audit.tsv")

def missing(series):
    return series.fillna("").astype(str).str.strip().str.lower().isin(
        {"", "na", "n/a", "none", "not provided", "-"}
    )

def main():
    if not INPUT.exists():
        raise FileNotFoundError(f"Input not found: {INPUT}")

    df = pd.read_csv(INPUT, sep="\t", dtype=str, keep_default_na=False)

    required = [
        "GeneSymbol", "VariationID", "Name", "ClinicalSignificance",
        "ReviewStatus", "NumberSubmitters", "LastEvaluated", "Assembly",
        "Chromosome", "PositionVCF", "ReferenceAlleleVCF",
        "AlternateAlleleVCF", "RCVaccession", "PhenotypeList",
        "SCVsForAggregateGermlineClassification",
    ]
    absent = sorted(set(required) - set(df.columns))
    if absent:
        raise ValueError(f"Missing columns: {absent}")

    df["SourceClassification"] = df["ClinicalSignificance"]

    df["Flag_ExplicitConflict"] = df["ClinicalSignificance"].str.contains(
        "conflicting", case=False, na=False
    )
    df["Flag_LimitedReviewStatus"] = df["ReviewStatus"].str.contains(
        "no assertion criteria provided|no assertion provided",
        case=False, regex=True, na=False
    )
    df["Flag_NoSubmitterCount"] = (
        missing(df["NumberSubmitters"])
        | ~df["NumberSubmitters"].str.fullmatch(r"\d+")
    )
    df["Flag_NoEvaluationDate"] = missing(df["LastEvaluated"])
    df["Flag_NoPhenotype"] = missing(df["PhenotypeList"])
    df["Flag_NoSCVReferences"] = missing(
        df["SCVsForAggregateGermlineClassification"]
    )
    df["Flag_NoVCFAlleles"] = (
        missing(df["ReferenceAlleleVCF"])
        | missing(df["AlternateAlleleVCF"])
    )

    flags = [
        "Flag_ExplicitConflict", "Flag_LimitedReviewStatus",
        "Flag_NoSubmitterCount", "Flag_NoEvaluationDate",
        "Flag_NoPhenotype", "Flag_NoSCVReferences", "Flag_NoVCFAlleles",
    ]
    df["AuditFlagCount"] = df[flags].astype(bool).sum(axis=1)

    df["ReviewPriority"] = "Routine review"
    df.loc[df["AuditFlagCount"] >= 2, "ReviewPriority"] = "Evidence review"
    df.loc[df["Flag_ExplicitConflict"], "ReviewPriority"] = (
        "Conflict review required"
    )

    columns = [
        "GeneSymbol", "VariationID", "Name", "Assembly", "Chromosome",
        "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF",
        "SourceClassification", "ReviewStatus", "NumberSubmitters",
        "LastEvaluated", "RCVaccession", "PhenotypeList",
        "SCVsForAggregateGermlineClassification", *flags,
        "AuditFlagCount", "ReviewPriority",
    ]

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    df[columns].to_csv(OUTPUT, sep="\t", index=False)

    print("EvidenceTrace audit created.")
    print(f"Rows audited: {len(df):,}")
    print(f"Output: {OUTPUT}")
    print("\nReview priority:")
    print(df["ReviewPriority"].value_counts().to_string())
    print("\nAudit flag counts:")
    print(df[flags].astype(bool).sum().sort_values(ascending=False).to_string())

if __name__ == "__main__":
    main()
