
import pytest

from geneweaver.gpu.alignment import encode_batch


def test_encode_batch():
    pairs = [
        ("ACGT", "TGCA"),
        ("AAAA", "CCCC")
    ]

    a, b, len_a, len_b = encode_batch(pairs)

    assert a.shape == (2, 256)
    assert b.shape == (2, 256)

    assert a[0, :4].tolist() == [0, 1, 2, 3]
    assert b[0, :4].tolist() == [3, 2, 1, 0]

    assert len_a.tolist() == [4, 4]
    assert len_b.tolist() == [4, 4]


def test_encode_ambiguous_base():
    a, b, len_a, len_b = encode_batch(
        [("ACNT", "ACGT")]
    )

    assert a[0, :4].tolist() == [0, 1, 4, 3]


def test_reject_long_sequence():
    with pytest.raises(ValueError):
        encode_batch([
            ("A" * 257, "ACGT")
        ])


def test_reject_empty_sequence():
    with pytest.raises(ValueError):
        encode_batch([
            ("", "ACGT")
        ])


def test_reject_invalid_base():
    with pytest.raises(ValueError):
        encode_batch([
            ("ACGX", "ACGT")
        ])
