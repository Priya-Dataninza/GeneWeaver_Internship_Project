from Bio import SeqIO
import numpy as np

def parse_fasta(file_path):
    """
    Read a FASTA file using BioPython.

    Returns:
        dict: {sequence_header: DNA_sequence}
    """
    sequences = {}

    for record in SeqIO.parse(file_path, "fasta"):
        header = record.description
        sequence = str(record.seq).upper()

        if sequence:
            sequences[header] = sequence

    if not sequences:
        raise ValueError("No valid FASTA sequences found.")

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

DNA_ENCODING = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
    "N": 4
}


def encode_dna(sequence):
    """Convert DNA sequence into a NumPy uint8 array."""

    sequence = sequence.upper()

    invalid = set(sequence) - set(DNA_ENCODING)

    if invalid:
        raise ValueError(f"Invalid DNA characters: {invalid}")

    return np.fromiter(
        (DNA_ENCODING[base] for base in sequence),
        dtype=np.uint8,
        count=len(sequence)
    )


def parse_fasta_array_chunks(file_path, chunk_size=100_000):
    """Generate encoded NumPy arrays from FASTA chunks."""

    for chunk in parse_fasta_in_chunks(file_path, chunk_size):
        yield {
            "header": chunk["header"],
            "start": chunk["start"],
            "sequence": encode_dna(chunk["sequence"])
        }