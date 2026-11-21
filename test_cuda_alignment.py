
import pytest

from geneweaver.cpu.alignment import needleman_wunsch
from geneweaver.gpu.alignment import run_cuda_alignment


# CUDA is optional on the development laptop.
try:
    from numba import cuda
    CUDA_AVAILABLE = cuda.is_available()
except Exception:
    CUDA_AVAILABLE = False


@pytest.mark.skipif(
    not CUDA_AVAILABLE,
    reason="CUDA GPU is not available"
)
def test_cuda_alignment_matches_cpu():
    """Verify that CUDA and CPU alignment scores match."""

    sequence_pairs = [
        ("ACGT", "ACGT"),
        ("AAAA", "AAAT"),
        ("ACGTACGT", "TGCATGCA"),
        ("GGGGCCCC", "GGGGAAAA"),
        ("ATCGATCG", "ATCGATCG"),
    ]

    # Calculate reference scores using CPU.
    cpu_scores = [
        needleman_wunsch(seq_a, seq_b)
        for seq_a, seq_b in sequence_pairs
    ]

    # Calculate scores using CUDA.
    gpu_scores = run_cuda_alignment(sequence_pairs)

    # Both implementations must produce identical results.
    assert gpu_scores == cpu_scores


@pytest.mark.skipif(
    not CUDA_AVAILABLE,
    reason="CUDA GPU is not available"
)
def test_cuda_alignment_batch():
    """Verify CUDA alignment for a larger batch."""

    sequence_pairs = [
        ("ACGT" * 16, "TGCA" * 16)
        for _ in range(100)
    ]

    cpu_scores = [
        needleman_wunsch(a, b)
        for a, b in sequence_pairs
    ]

    gpu_scores = run_cuda_alignment(sequence_pairs)

    assert len(gpu_scores) == 100
    assert gpu_scores == cpu_scores
