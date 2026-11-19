
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



