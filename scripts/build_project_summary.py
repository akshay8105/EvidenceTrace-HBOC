from pathlib import Path
import pandas as pd

ROOT = Path(".")
OUT = ROOT / "results/reports/project_quality_summary.tsv"

def main():
    audit_path = ROOT / "results/reports/evidencetrace_evidence_audit.tsv"
    candidates_path = ROOT / "results/reports/refined_disagreement_candidates.tsv"
    provenance_path = ROOT / "results/reports/evidence_provenance_report.tsv"

    for path in (audit_path, candidates_path, provenance_path):
        if not path.exists():
            raise FileNotFoundError(f"Missing required report: {path}")

    audit = pd.read_csv(audit_path, sep="\t", dtype=str, keep_default_na=False)
    candidates = pd.read_csv(
        candidates_path, sep="\t", dtype=str, keep_default_na=False
    )
    provenance = pd.read_csv(
        provenance_path, sep="\t", dtype=str, keep_default_na=False
    )

    metrics = [
        ("aggregate_audit_records", len(audit),
         "Rows in the rule-based aggregate evidence audit"),
        ("refined_disagreement_candidate_groups", len(candidates),
         "Variant-condition groups flagged for review"),
        ("provenance_report_rows", len(provenance),
         "Rows in the provenance report"),
        ("provenance_report_columns", len(provenance.columns),
         "Columns in the provenance report"),
        ("duplicate_variant_condition_rows",
         int(provenance.duplicated(
             subset=["VariationID", "ReportedPhenotypeInfo"]
         ).sum()),
         "Duplicate variant-condition rows in provenance report"),
    ]

    for priority, count in provenance["ReviewPriority"].value_counts().items():
        key = priority.lower().replace(" ", "_")
        metrics.append((f"priority_{key}", int(count),
                        "Rule-based review-priority group count"))

    summary = pd.DataFrame(
        metrics, columns=["Metric", "Value", "Interpretation"]
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(OUT, sep="\t", index=False)

    print("Project quality summary created.")
    print(summary.to_string(index=False))
    print(f"\nSaved to: {OUT}")

if __name__ == "__main__":
    main()
