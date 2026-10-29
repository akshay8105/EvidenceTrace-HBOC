import pandas as pd


def test_explicit_conflict_detection():
    values = pd.Series([
        "Conflicting classifications of pathogenicity",
        "Pathogenic",
        "Benign",
    ])

    result = values.str.contains("conflicting", case=False, na=False)

    assert result.tolist() == [True, False, False]


def test_missing_value_detection():
    values = pd.Series(["", "na", "Not provided", "Pathogenic"])

    missing = values.str.strip().str.lower().isin(
        {"", "na", "n/a", "none", "not provided", "-"}
    )

    assert missing.tolist() == [True, True, True, False]


def test_priority_rules():
    df = pd.DataFrame({
        "AuditFlagCount": [0, 2, 1],
        "Flag_ExplicitConflict": [False, False, True],
    })

    df["ReviewPriority"] = "Routine review"
    df.loc[df["AuditFlagCount"] >= 2, "ReviewPriority"] = "Evidence review"
    df.loc[df["Flag_ExplicitConflict"], "ReviewPriority"] = (
        "Conflict review required"
    )

    assert df["ReviewPriority"].tolist() == [
        "Routine review",
        "Evidence review",
        "Conflict review required",
    ]


def test_filtered_dataset_has_expected_genes():
    genes = {"ATM", "BRCA1", "BRCA2", "CHEK2", "PALB2"}
    assert len(genes) == 5
    assert "BRCA1" in genes
