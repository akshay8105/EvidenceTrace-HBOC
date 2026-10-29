# EvidenceTrace-HBOC

A reproducible bioinformatics workflow for auditing ClinVar variant interpretations and tracing submission-level evidence in hereditary breast and ovarian cancer susceptibility genes.

## Project overview

This research prototype analyzes ClinVar aggregate records and individual submissions for ATM, BRCA1, BRCA2, CHEK2, and PALB2.

It audits clinical-significance labels, review statuses, and submission metadata. It identifies variant-condition groups with differing normalized classifications and creates reports to support manual evidence review.

## Current results

| Metric | Result |
|---|---:|
| Aggregate audit records | 68,198 |
| Disagreement candidate groups | 10,738 |
| Provenance report rows | 10,738 |
| Provenance report columns | 22 |
| Duplicate variant-condition rows | 0 |

These are pipeline output counts, not necessarily unique variants. Review priorities are heuristic flags, not clinical risk scores or validated pathogenicity predictions.

## Workflow

1. Filter ClinVar records for five selected genes and GRCh38.
2. Audit aggregate significance labels and review statuses.
3. Extract individual submissions using VariationID.
4. Normalize classification labels while retaining original labels.
5. Identify candidate disagreements by variant and reported condition.
6. Assign rule-based review priorities.
7. Attach submission provenance, including SCVs, submitters, review statuses, and evaluation dates.
8. Validate report structure and summarize results.

## Key outputs

- `results/reports/evidencetrace_evidence_audit.tsv`
- `results/reports/refined_disagreement_candidates.tsv`
- `results/reports/prioritized_disagreement_review.tsv`
- `results/reports/evidence_provenance_report.tsv`
- `results/reports/project_quality_summary.tsv`

Audit plots are stored in `results/plots/`.

## Limitations

- ClinVar interpretations may change over time.
- Different submissions may use different conditions and evidence.
- Submission counts do not represent independent evidence sources.
- Disagreement flags require manual review.
- The priority rules have not been clinically validated.
- This prototype does not independently establish ACMG/AMP classifications or clinical utility.

## Reproducibility

Use Python 3.12 or a compatible version. Activate the virtual environment with `source .venv/bin/activate`, then run the workflow scripts in dependency order and validate with `pytest -q`.

## Disclaimer

For educational and research use only. Not intended for clinical decision-making, diagnosis, or treatment selection.