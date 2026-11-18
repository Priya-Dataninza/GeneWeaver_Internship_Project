
from geneweaver.cpu.baseline import (
    reverse_complement,
    search_both_strands,
)


def test_reverse_complement():
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AAGC") == "GCTT"


def test_forward_strand_search():
    results = search_both_strands(
        "ACGTACGT",
        "TTTACGTACGTAGGTTT",
        max_mismatches=2
    )

    assert any(result["strand"] == "+" for result in results)


def test_reverse_strand_search():
    results = search_both_strands(
        "ACGTACGT",
        "CCTACGTACGT",
        max_mismatches=2
    )

    assert any(
        result["strand"] == "-"
        and result["position"] == 0
        for result in results
    )
