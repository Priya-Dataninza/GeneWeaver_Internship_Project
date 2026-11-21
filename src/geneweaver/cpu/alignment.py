
import time


def needleman_wunsch(
    sequence_a,
    sequence_b,
    match=2,
    mismatch=-1,
    gap=-2
):
    """
    Global DNA sequence alignment using Needleman-Wunsch.
    Returns the optimal alignment score.
    """

    sequence_a = sequence_a.upper()
    sequence_b = sequence_b.upper()

    if not sequence_a or not sequence_b:
        raise ValueError("Sequences cannot be empty.")

    valid_bases = {"A", "C", "G", "T", "N"}

    if set(sequence_a) - valid_bases:
        raise ValueError("Invalid DNA bases in sequence A.")

    if set(sequence_b) - valid_bases:
        raise ValueError("Invalid DNA bases in sequence B.")

    rows = len(sequence_a) + 1
    cols = len(sequence_b) + 1

    matrix = [
        [0] * cols
        for _ in range(rows)
    ]

    # Initialize first column
    for i in range(rows):
        matrix[i][0] = i * gap

    # Initialize first row
    for j in range(cols):
        matrix[0][j] = j * gap

    # Fill dynamic programming matrix
    for i in range(1, rows):
        for j in range(1, cols):

            if sequence_a[i - 1] == sequence_b[j - 1]:
                score = match
            else:
                score = mismatch

            diagonal = matrix[i - 1][j - 1] + score
            up = matrix[i - 1][j] + gap
            left = matrix[i][j - 1] + gap

            matrix[i][j] = max(
                diagonal,
                up,
                left
            )

    return matrix[-1][-1]


def benchmark_alignment(sequence_a, sequence_b):
    """Measure CPU alignment execution time."""

    start_time = time.perf_counter()

    score = needleman_wunsch(
        sequence_a,
        sequence_b
    )

    elapsed_time = time.perf_counter() - start_time

    return {
        "alignment_score": score,
        "execution_time_seconds": elapsed_time
    }
