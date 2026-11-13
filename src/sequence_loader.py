def load_sequences(file_path):
    """Read DNA sequences from a text file."""

    sequences = []

    try:
        with open(file_path, "r") as file:
            for line in file:
                sequence = line.strip().upper()

                if sequence:
                    sequences.append(sequence)

        return sequences

    except FileNotFoundError:
        print("Error: Input file not found")
        return []
        print("Total sequences loaded:", len(sequences))

        for sequence in sequences:
            print("Sequence:", sequence)
            print("Length:", len(sequence))