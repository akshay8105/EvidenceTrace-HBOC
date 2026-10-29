from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

REPORTS = Path("results/reports")
PLOTS = Path("results/plots")
PLOTS.mkdir(parents=True, exist_ok=True)

audit_path = REPORTS / "evidencetrace_evidence_audit.tsv"
if not audit_path.exists():
    raise FileNotFoundError(f"Run the evidence audit first: {audit_path}")

df = pd.read_csv(audit_path, sep="\t", dtype=str)

# Convert flags explicitly from text to Boolean values.
df["Flag_ExplicitConflict"] = (
    df["Flag_ExplicitConflict"].str.lower().eq("true")
)

# Chart 1: review-priority counts by gene.
counts = pd.crosstab(df["GeneSymbol"], df["ReviewPriority"])
priority_order = [
    "Routine review",
    "Evidence review",
    "Conflict review required",
]
counts = counts.reindex(columns=priority_order, fill_value=0)

ax = counts.plot(kind="bar", stacked=True, figsize=(10, 6))
ax.set_title("EvidenceTrace-HBOC: Audit Review Categories")
ax.set_xlabel("Gene")
ax.set_ylabel("ClinVar summary rows")
ax.legend(title="Rule-based review category")
plt.xticks(rotation=0)
plt.tight_layout()
plt.savefig(PLOTS / "review_categories_by_gene.png", dpi=200)
plt.close()

# Chart 2: most common source clinical-significance labels.
top_labels = df["SourceClassification"].value_counts().head(12).sort_values()

fig, ax = plt.subplots(figsize=(11, 7))
top_labels.plot(kind="barh", ax=ax)
ax.set_title("ClinVar Source Classification Labels")
ax.set_xlabel("Number of rows")
ax.set_ylabel("Source classification")
plt.tight_layout()
plt.savefig(PLOTS / "source_classifications.png", dpi=200)
plt.close()

# Chart 3: explicit conflict-flag percentage by gene.
conflict_rate = (
    df.groupby("GeneSymbol")["Flag_ExplicitConflict"].mean() * 100
).sort_values()

fig, ax = plt.subplots(figsize=(9, 5))
conflict_rate.plot(kind="bar", ax=ax)
ax.set_title("Rows Explicitly Labelled as Conflicting by Gene")
ax.set_xlabel("Gene")
ax.set_ylabel("Rows with explicit conflict label (%)")
plt.xticks(rotation=0)
plt.tight_layout()
plt.savefig(PLOTS / "explicit_conflict_percentage.png", dpi=200)
plt.close()

print("Charts created:")
for path in sorted(PLOTS.glob("*.png")):
    print(f"{path} ({path.stat().st_size:,} bytes)")
