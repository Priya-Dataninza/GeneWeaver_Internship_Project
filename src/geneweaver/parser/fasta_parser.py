
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


def chunk_sequence(sequence, chunk_size=100_000):
    """
    Divide a DNA sequence into smaller chunks.

    Args:
        sequence: DNA sequence string.
        chunk_size: Maximum bases per chunk.

    Returns:
        List of DNA sequence chunks.
    """
    if chunk_size <= 0:
        raise ValueError("Chunk size must be positive.")

    sequence = sequence.upper()

    for start in range(0, len(sequence), chunk_size):
        yield {
            "start": start,
            "sequence": sequence[start:start + chunk_size]
        }


def parse_fasta_in_chunks(file_path, chunk_size=100_000):
    """
    Parse a FASTA file and yield its sequences in chunks.
    """
    sequences = parse_fasta(file_path)

    for header, sequence in sequences.items():
        for chunk in chunk_sequence(sequence, chunk_size):
            yield {
                "header": header,
                "start": chunk["start"],
                "sequence": chunk["sequence"]
            }

