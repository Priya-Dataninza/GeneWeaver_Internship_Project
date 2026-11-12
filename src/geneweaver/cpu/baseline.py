
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
