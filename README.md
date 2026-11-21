"# GeneWeaver_Internship_Project" 
"# GeneWeaver_Internship_Project" 

## CUDA Acceleration and Performance Benchmark

GeneWeaver uses Numba CUDA to accelerate DNA sequence alignment through GPU parallelization.

### Implementation
- Encodes DNA bases into numerical representations.
- Uses a CUDA kernel to calculate Needleman–Wunsch alignment scores.
- Processes multiple sequence pairs in parallel.
- Uses rolling dynamic programming arrays to reduce memory requirements.
- Validates GPU results against the CPU implementation.

### Initial Benchmark Results

| Metric | Result |
|---|---:|
| Sequence pairs | 1,000 |
| Sequence length | 64 bases |
| CPU execution time | 1.842064 seconds |
| GPU execution time | 0.282563 seconds |
| Measured speedup | 6.52× |
| CPU/GPU correctness | Passed |

### Performance Goal

The project targets a 100× GPU speedup. The current measured result is 6.52× on the initial benchmark workload.

Further optimization will focus on:
- Increasing batch size and GPU utilization.
- Reducing host-to-device and device-to-host transfer overhead.
- Optimizing CUDA kernel execution.
- Evaluating performance across larger workloads.

**Note:** The 6.52× result is an initial benchmark, not proof that the 100× target has been achieved.
