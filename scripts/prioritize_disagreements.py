from pathlib import Path
import pandas as pd

INPUT = Path("results/reports/refined_disagreement_candidates.tsv")
OUTPUT = Path("results/reports/prioritized_disagreement_review.tsv")

def assign_priority(value):
    labels = set(str(value).split(" | "))

    pathogenic = labels.intersection({
        "Pathogenic", "Likely pathogenic",
        "Pathogenic/Likely pathogenic"
    })
    benign = labels.intersection({
        "Benign", "Likely benign",
        "Benign/Likely benign"
    })
    uncertain = "Uncertain significance" in labels

    if pathogenic and benign:
        return "High review priority"
    if uncertain and (pathogenic or benign):
        return "Moderate review priority"
    return "Routine follow-up"

def main():
    if not INPUT.exists():
        raise FileNotFoundError(f"Missing input report: {INPUT}")

    df = pd.read_csv(INPUT, sep="\t", dtype=str, keep_default_na=False)

    required = {
        "VariationID", "MappedGeneSymbol",
        "NormalizedClassifications", "ReportedPhenotypeInfo"
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")

    df["ReviewPriority"] = df["NormalizedClassifications"].apply(
        assign_priority
    )

    priority_order = {
        "High review priority": 0,
        "Moderate review priority": 1,
        "Routine follow-up": 2
    }
    df["_order"] = df["ReviewPriority"].map(priority_order)
    df = df.sort_values(
        ["_order", "MappedGeneSymbol", "VariationID"]
    ).drop(columns="_order")

    df.to_csv(OUTPUT, sep="\t", index=False)

    print("Prioritized review report created.")
    print(f"Candidate groups: {len(df):,}")
    print("\nPriority counts:")
    print(df["ReviewPriority"].value_counts().to_string())
    print(f"\nOutput: {OUTPUT}")

if __name__ == "__main__":
    main()
