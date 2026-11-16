
import numpy as np

# DNA encoding:
# A = 0, C = 1, G = 2, T = 3
DNA_ENCODING = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3
}


def encode_dna(sequence):
    """Convert DNA letters into numerical values."""
    sequence = sequence.upper()

    try:
        return np.array(
            [DNA_ENCODING[base] for base in sequence],
            dtype=np.uint8
        )
    except KeyError as error:
        raise ValueError(
            f"Invalid DNA character: {error.args[0]}"
        ) from error


def create_cuda_kernel():
    """Create the Numba CUDA kernel when CUDA is available."""
    try:
        from numba import cuda
    except ImportError as error:
        raise RuntimeError(
            "Numba CUDA could not be imported. "
            "Use a compatible GPU environment such as Colab."
        ) from error

    @cuda.jit
    def off_target_kernel(
        guide,
        candidates,
        mismatch_counts,
        valid_matches,
        guide_length,
        max_mismatches
    ):
        index = cuda.grid(1)

        if index >= candidates.shape[0]:
            return

        mismatches = 0

        # Compare guide and candidate bases
        for j in range(guide_length):
            if guide[j] != candidates[index, j]:
                mismatches += 1

        # NGG PAM: N can be any base; last two must be G
        pam_valid = (
            candidates[index, guide_length + 1] == 2
            and candidates[index, guide_length + 2] == 2
        )

        mismatch_counts[index] = mismatches

        if pam_valid and mismatches <= max_mismatches:
            valid_matches[index] = 1
        else:
            valid_matches[index] = 0

    return off_target_kernel


def run_gpu_search(guide, candidates, max_mismatches=2):
    """
    Transfer input arrays to GPU, execute the CUDA kernel,
    and transfer the results back to the host.
    """
    from numba import cuda

    guide_array = encode_dna(guide)

    candidate_array = np.array(
        [encode_dna(candidate) for candidate in candidates],
        dtype=np.uint8
    )

    if candidate_array.ndim != 2:
        raise ValueError("Candidates must have equal lengths.")

    if candidate_array.shape[1] != len(guide) + 3:
        raise ValueError(
            "Each candidate must contain guide-length DNA "
            "followed by a 3-base PAM."
        )

    if len(candidate_array) == 0:
        return [], []

    # Host-to-device transfers
    device_guide = cuda.to_device(guide_array)
    device_candidates = cuda.to_device(candidate_array)

    mismatch_counts = cuda.device_array(
        len(candidates), dtype=np.int32
    )
    valid_matches = cuda.device_array(
        len(candidates), dtype=np.int32
    )

    kernel = create_cuda_kernel()

    threads_per_block = 128
    blocks_per_grid = (
        len(candidates) + threads_per_block - 1
    ) // threads_per_block

    kernel[blocks_per_grid, threads_per_block](
        device_guide,
        device_candidates,
        mismatch_counts,
        valid_matches,
        len(guide),
        max_mismatches
    )

    # Device-to-host transfers
    host_mismatches = mismatch_counts.copy_to_host()
    host_matches = valid_matches.copy_to_host()

    return host_mismatches, host_matches
