import time


def needleman_wunsch(
    sequence_a,
    sequence_b,
    match=2,
    mismatch=-1,
    gap=-2
):
    """
    Perform global DNA sequence alignment
    using the Needleman-Wunsch algorithm.
    """

    sequence_a = validate_dna(sequence_a)
    sequence_b = validate_dna(sequence_b)

    rows = len(sequence_a) + 1
    cols = len(sequence_b) + 1

    # Initialize scoring matrix
    matrix = [
        [0] * cols
        for _ in range(rows)
    ]

    for i in range(rows):
        matrix[i][0] = i * gap

    for j in range(cols):
        matrix[0][j] = j * gap

    # Dynamic programming
    for i in range(1, rows):
        for j in range(1, cols):

            score = (
                match
                if sequence_a[i - 1] == sequence_b[j - 1]
                else mismatch
            )

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
    """Measure Needleman-Wunsch execution time."""

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


def count_mismatches(guide, candidate):
    """Count differences between guide and candidate DNA."""
    return sum(
        base1 != base2
        for base1, base2 in zip(guide, candidate)
    )


def has_ngg_pam(pam):
    """Check the SpCas9 NGG PAM pattern."""
    return (
        len(pam) == 3
        and pam[1] == "G"
        and pam[2] == "G"
    )


def search_off_targets(
    guide,
    genome,
    max_mismatches=2
):
    """
    Search a DNA sequence for guide-length candidates
    followed by an NGG PAM.

    Returns candidate positions and mismatch counts.
    """
    guide = guide.upper()
    genome = genome.upper()

    guide_length = len(guide)
    pam_length = 3
    window_length = guide_length + pam_length

    results = []

    for position in range(len(genome) - window_length + 1):
        candidate = genome[position:position + guide_length]
        pam = genome[
            position + guide_length:
            position + window_length
        ]

        if not has_ngg_pam(pam):
            continue

        mismatches = count_mismatches(guide, candidate)

        if mismatches <= max_mismatches:
            results.append({
                "position": position,
                "candidate": candidate,
                "pam": pam,
                "mismatches": mismatches
            })

    return results



def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence."""
    complement = str.maketrans(
        "ACGTN",
        "TGCAN"
    )
    return sequence.upper().translate(complement)[::-1]


def validate_dna(sequence):
    """Validate that a DNA sequence contains only A, C, G, T."""
    sequence = sequence.upper()

    if not sequence:
        raise ValueError("DNA sequence cannot be empty.")

    invalid_bases = set(sequence) - {"A", "C", "G", "T"}

    if invalid_bases:
        raise ValueError(
            f"Invalid DNA bases found: {sorted(invalid_bases)}"
        )

    return sequence

def search_both_strands(
    guide,
    genome,
    max_mismatches=2
):
    """Search both forward and reverse DNA strands."""

    guide = validate_dna(guide)
    genome = validate_dna(genome)
    # Search the forward strand
    forward_results = search_off_targets(
        guide,
        genome,
        max_mismatches
    )

    for result in forward_results:
        result["strand"] = "+"

    # Search the reverse-complement strand
    reverse_genome = reverse_complement(genome)

    reverse_results = search_off_targets(
        guide,
        reverse_genome,
        max_mismatches
    )

    # Convert reverse-oriented positions to reference coordinates
    window_length = len(guide) + 3

    for result in reverse_results:
        result["position"] = (
            len(genome)
            - result["position"]
            - window_length
        )
        result["strand"] = "-"

    return sorted(
        forward_results + reverse_results,
        key=lambda result: result["position"]
    )



