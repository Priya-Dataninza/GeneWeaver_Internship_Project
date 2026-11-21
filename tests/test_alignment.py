
import pytest

from geneweaver.cpu.alignment import needleman_wunsch


def test_identical_sequences():
    assert needleman_wunsch("ACGT", "ACGT") == 8


def test_alignment_with_gap():
    assert needleman_wunsch("ACGT", "AGT") == 4


def test_alignment_with_mismatch():
    assert needleman_wunsch("ACGT", "AGGT") == 5


def test_empty_sequence():
    with pytest.raises(ValueError):
        needleman_wunsch("", "ACGT")


def test_invalid_dna():
    with pytest.raises(ValueError):
        needleman_wunsch("ACGX", "ACGT")
