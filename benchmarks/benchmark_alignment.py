from geneweaver.cpu.baseline import benchmark_alignment

sequence_a = "ACGTACGTACGTACGT"
sequence_b = "ACGTTCGTACGTACGT"

result = benchmark_alignment(sequence_a, sequence_b)

print("GENEWEAVER CPU ALIGNMENT")
print("------------------------")
print("Alignment score:", result["alignment_score"])
print("Execution time:", result["execution_time_seconds"], "seconds")