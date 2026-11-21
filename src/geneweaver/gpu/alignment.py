
import numpy as np

DNA_ENCODING = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
    "N": 4
}

MAX_SEQ_LEN = 256


def encode_batch(sequence_pairs):
    """Encode DNA sequence pairs into fixed-width NumPy arrays."""

    count = len(sequence_pairs)

    sequences_a = np.zeros(
        (count, MAX_SEQ_LEN), dtype=np.uint8
    )

    sequences_b = np.zeros(
        (count, MAX_SEQ_LEN), dtype=np.uint8
    )

    lengths_a = np.zeros(count, dtype=np.int32)
    lengths_b = np.zeros(count, dtype=np.int32)

    for index, (seq_a, seq_b) in enumerate(sequence_pairs):

        seq_a = seq_a.upper()
        seq_b = seq_b.upper()

        if not seq_a or not seq_b:
            raise ValueError("Sequences cannot be empty.")

        if len(seq_a) > MAX_SEQ_LEN or len(seq_b) > MAX_SEQ_LEN:
            raise ValueError(
                f"Sequences must not exceed {MAX_SEQ_LEN} bases."
            )

        for seq in (seq_a, seq_b):
            invalid = set(seq) - set(DNA_ENCODING)

            if invalid:
                raise ValueError(
                    f"Invalid DNA bases: {sorted(invalid)}"
                )

        sequences_a[index, :len(seq_a)] = [
            DNA_ENCODING[base] for base in seq_a
        ]

        sequences_b[index, :len(seq_b)] = [
            DNA_ENCODING[base] for base in seq_b
        ]

        lengths_a[index] = len(seq_a)
        lengths_b[index] = len(seq_b)

    return sequences_a, sequences_b, lengths_a, lengths_b


def create_alignment_kernel():
    """Create the Numba CUDA alignment kernel lazily."""

    from numba import cuda, int32

    @cuda.jit
    def alignment_kernel(
        sequences_a,
        sequences_b,
        lengths_a,
        lengths_b,
        output,
        match_score,
        mismatch_score,
        gap_penalty
    ):

        index = cuda.grid(1)

        if index >= sequences_a.shape[0]:
            return

        len_a = lengths_a[index]
        len_b = lengths_b[index]

        # Rolling dynamic programming row
        dp = cuda.local.array(
            257,
            dtype=int32
        )

        for j in range(len_b + 1):
            dp[j] = j * gap_penalty

        for i in range(1, len_a + 1):

            previous_diagonal = dp[0]
            dp[0] = i * gap_penalty

            for j in range(1, len_b + 1):

                old_value = dp[j]

                if (
                    sequences_a[index, i - 1]
                    == sequences_b[index, j - 1]
                ):
                    score = match_score
                else:
                    score = mismatch_score

                diagonal = previous_diagonal + score
                up = old_value + gap_penalty
                left = dp[j - 1] + gap_penalty

                best = diagonal

                if up > best:
                    best = up

                if left > best:
                    best = left

                dp[j] = best
                previous_diagonal = old_value

        output[index] = dp[len_b]

    return alignment_kernel


def run_cuda_alignment(
    sequence_pairs,
    match_score=2,
    mismatch_score=-1,
    gap_penalty=-2
):
    """Transfer DNA sequences to GPU and calculate alignment scores."""

    if not sequence_pairs:
        return []

    if match_score <= 0:
        raise ValueError("Match score must be positive.")

    if mismatch_score > 0 or gap_penalty > 0:
        raise ValueError(
            "Mismatch score and gap penalty must be non-positive."
        )

    from numba import cuda

    if not cuda.is_available():
        raise RuntimeError(
            "CUDA GPU is unavailable. Use an NVIDIA GPU environment."
        )

    host_a, host_b, lengths_a, lengths_b = encode_batch(
        sequence_pairs
    )

    # Host-to-device transfers
    device_a = cuda.to_device(host_a)
    device_b = cuda.to_device(host_b)
    device_lengths_a = cuda.to_device(lengths_a)
    device_lengths_b = cuda.to_device(lengths_b)

    device_output = cuda.device_array(
        len(sequence_pairs),
        dtype=np.int32
    )

    kernel = create_alignment_kernel()

    threads_per_block = 128

    blocks_per_grid = (
        len(sequence_pairs) + threads_per_block - 1
    ) // threads_per_block

    kernel[blocks_per_grid, threads_per_block](
        device_a,
        device_b,
        device_lengths_a,
        device_lengths_b,
        device_output,
        match_score,
        mismatch_score,
        gap_penalty
    )

    # Device-to-host transfer
    host_output = device_output.copy_to_host()

    return host_output.tolist()
