import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from refine_disagreements import normalize_classification


def test_pathogenic_labels():
    assert normalize_classification("Pathogenic") == "Pathogenic"
    assert normalize_classification("Likely pathogenic") == "Likely pathogenic"
    assert normalize_classification("Pathogenic/Likely pathogenic") == (
        "Pathogenic/Likely pathogenic"
    )


def test_benign_labels():
    assert normalize_classification("Benign") == "Benign"
    assert normalize_classification("Likely benign") == "Likely benign"
    assert normalize_classification("Benign/Likely benign") == (
        "Benign/Likely benign"
    )


def test_uncertain_significance():
    assert normalize_classification("Uncertain significance") == (
        "Uncertain significance"
    )


def test_whitespace_and_case():
    assert normalize_classification("  LIKELY BENIGN  ") == "Likely benign"


def test_missing_labels():
    assert normalize_classification("") == "Missing/unspecified"
    assert normalize_classification("-") == "Missing/unspecified"
    assert normalize_classification("not provided") == "Missing/unspecified"


def test_unmapped_label_is_preserved():
    assert normalize_classification("Some new classification") == (
        "Other/unmapped: some new classification"
    )
