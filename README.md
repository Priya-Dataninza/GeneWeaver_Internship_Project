# GeneWeaver – GPU-Accelerated DNA Sequence Alignment

GeneWeaver is a Bioinformatics and High-Performance Computing (HPC) project designed to accelerate DNA sequence analysis using CPU and GPU computing.

The project explores parallel computing techniques to improve the performance of DNA sequence alignment, with a long-term goal of accelerating CRISPR off-target sequence analysis.

## Project Objectives

- Parse DNA sequences from FASTA files.
- Divide large sequences into manageable chunks.
- Implement a CPU-based sequence alignment baseline.
- Accelerate alignment using NVIDIA CUDA and Numba.
- Compare CPU and GPU execution times.
- Validate GPU results against CPU results.
- Work towards a target of 100× GPU acceleration.

## Technology Stack

| Technology | Purpose |
|---|---|
| Python | Core implementation |
| NumPy | Numerical arrays and DNA encoding |
| Numba CUDA | GPU kernel development |
| BioPython | Biological sequence processing |
| FASTA | DNA sequence input format |
| Pytest | Automated testing |
| Google Colab | NVIDIA GPU execution and benchmarking |

## Project Structure

```text
Geneweaver/
│
├── benchmarks/
│   ├── benchmark_alignment.py
│   ├── benchmark_cuda_alignment.py
│   └── __init__.py
│
├── data/
│   └── sample_genome.fasta
│
├── src/
│   └── geneweaver/
│       ├── cpu/
│       │   ├── alignment.py
│       │   └── baseline.py
│       │
│       ├── gpu/
│       │   ├── alignment.py
│       │   └── kernels.py
│       │
│       ├── parser/
│       │   └── fasta_parser.py
│       │
│       ├── distributed/
│       │   └── scheduler.py
│       │
│       ├── ui/
│       └── main.py
│
└── tests/
    ├── test_alignment.py
    ├── test_cpu_baseline.py
    ├── test_fasta_chunks.py
    ├── test_gpu_encoding.py
    └── test_cuda_alignment.py
```

## Core Features

### 1. FASTA Parsing
- Reads DNA sequences from FASTA files.
- Supports sequence chunking.
- Prepares DNA data for downstream processing.

### 2. CPU Alignment
- Implements the Needleman–Wunsch global alignment algorithm.
- Uses dynamic programming.
- Provides a reference implementation for performance comparison.

### 3. GPU Acceleration
- Encodes DNA bases as numerical values:
  - A = 0
  - C = 1
  - G = 2
  - T = 3
  - N = 4
- Uses Numba CUDA to execute alignment calculations on an NVIDIA GPU.
- Uses rolling dynamic programming arrays to reduce per-thread memory requirements.
- Processes multiple sequence pairs in parallel.

### 4. Performance Benchmarking
- Measures CPU and GPU execution time.
- Compares results for correctness.
- Calculates GPU speedup.

## Initial CUDA Benchmark Results

The first successful CUDA benchmark processed 1,000 DNA sequence pairs, each containing 64 bases.

| Metric | Result |
|---|---:|
| Sequence pairs | 1,000 |
| Sequence length | 64 bases |
| CPU execution time | 1.842064 seconds |
| GPU execution time | 0.282563 seconds |
| Measured speedup | 6.52× |
| CPU/GPU result comparison | Passed |

**Initial benchmark speedup: 6.52×**

The project target is 100× GPU acceleration. This target has not yet been achieved.

Further optimization will investigate:
- Larger sequence batches.
- GPU thread and block utilization.
- Host-to-device and device-to-host transfer overhead.
- CUDA kernel optimization.
- Performance on larger DNA sequences.

*Benchmark timings depend on the GPU hardware, workload, and execution environment. The reported result is from the initial Google Colab benchmark.*

## Installation

Create and activate a Python virtual environment:

```bash
python -m venv .venv
```

Windows:

```bash
.venv\Scripts\activate
```

Install dependencies:

```bash
pip install numpy numba biopython pytest
```

For CUDA execution, use a compatible NVIDIA GPU environment with the required CUDA driver and Numba CUDA support. Google Colab with a GPU runtime can be used.

## Running Tests

Set the source directory in the Python path.

Windows PowerShell:

```powershell
$env:PYTHONPATH = "src"
```

Run:

```bash
python -m pytest tests/ -v
```

The CUDA correctness test is skipped automatically when CUDA is unavailable.

## Running the CUDA Benchmark

Run this in an NVIDIA GPU environment:

```bash
python benchmarks/benchmark_cuda_alignment.py
```

The benchmark compares CPU and GPU scores before reporting performance.

## Hardware Note

The development laptop uses Intel Iris Xe integrated graphics and does not provide an NVIDIA CUDA GPU. Therefore, CUDA execution and GPU benchmarking are performed using Google Colab's NVIDIA GPU runtime.

## Project Status

- [x] FASTA parsing and chunking
- [x] CPU alignment implementation
- [x] DNA numerical encoding
- [x] CUDA alignment kernel implementation
- [x] Initial CUDA execution
- [x] CPU/GPU correctness comparison
- [x] Initial performance benchmark
- [ ] Larger batch-size experiments
- [ ] CUDA performance optimization
- [ ] Achieve the 100× speedup target

---

**GeneWeaver – Exploring GPU Acceleration for Bioinformatics and High-Performance Computing.**