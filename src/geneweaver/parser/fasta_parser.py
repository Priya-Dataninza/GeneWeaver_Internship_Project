
def parse_fasta(file_path):
    """
    Read a FASTA file and return its DNA sequences.

    Returns:
        dict: {sequence_header: DNA_sequence}
    """
    sequences = {}
    current_header = None
    current_sequence = []

    with open(file_path, "r", encoding="utf-8") as fasta_file:
        for line in fasta_file:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_header is not None:
                    sequences[current_header] = "".join(current_sequence)

                current_header = line[1:].strip()
                current_sequence = []
            else:
                if current_header is None:
                    raise ValueError("Invalid FASTA: sequence found before header")

                current_sequence.append(line.upper())

    if current_header is not None:
        sequences[current_header] = "".join(current_sequence)

    return sequences
