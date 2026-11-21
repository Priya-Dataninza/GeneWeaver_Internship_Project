
import time

from geneweaver.cpu.alignment import needleman_wunsch
from geneweaver.gpu.alignment import run_cuda_alignment


def generate_pairs(count=1000):
    """Generate DNA sequence pairs for benchmarking."""

    base = "ACGT" * 16
    pairs = []

    for i in range(count):
        candidate = list(base)
        position = i % len(base)

        candidate[position] = (
            "A" if candidate[position] != "A" else "C"
        )

        pairs.append((base, "".join(candidate)))

    return pairs


def benchmark():
    pairs = generate_pairs(1000)

    # CPU benchmark
    start = time.perf_counter()

    cpu_scores = [
        needleman_wunsch(a, b)
        for a, b in pairs
    ]

    cpu_time = time.perf_counter() - start

    # GPU warm-up: includes first-time compilation
    run_cuda_alignment(pairs[:1])

    # GPU benchmark: includes transfers and execution
    start = time.perf_counter()

    gpu_scores = run_cuda_alignment(pairs)

    gpu_time = time.perf_counter() - start

    # Verify correctness
    assert cpu_scores == gpu_scores, "CPU and GPU scores differ!"

    speedup = cpu_time / gpu_time if gpu_time > 0 else 0

    print("\nGENEWEAVER CUDA BENCHMARK")
    print("------------------------")
    print("Sequence pairs:", len(pairs))
    print("CPU time:", round(cpu_time, 6), "seconds")
    print("GPU time:", round(gpu_time, 6), "seconds")
    print("Speedup:", round(speedup, 2), "x")
    print("CPU and GPU results match!")


if __name__ == "__main__":
    benchmark()

