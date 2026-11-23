import os
import sys
import tempfile
import time

from django.conf import settings
from django.shortcuts import render

# Add the src directory so Django can import GeneWeaver.
PROJECT_ROOT = settings.BASE_DIR
SRC_DIR = PROJECT_ROOT / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from geneweaver.parser.fasta_parser import parse_fasta
from geneweaver.cpu.baseline import search_both_strands


def home(request):
    context = {}

    if request.method == "POST":
        uploaded_file = request.FILES.get("fasta_file")
        guide = request.POST.get("guide", "").strip().upper()

        try:
            max_mismatches = int(
                request.POST.get("max_mismatches", 2)
            )

            if not uploaded_file:
                raise ValueError("Please upload a FASTA file.")

            if not uploaded_file.name.lower().endswith(
                (".fasta", ".fa", ".fna")
            ):
                raise ValueError("Please upload a valid FASTA file.")

            if not guide or set(guide) - {"A", "C", "G", "T"}:
                raise ValueError("Guide must contain only A, C, G and T.")

            if max_mismatches not in range(5):
                raise ValueError("Invalid mismatch limit.")

            # Save the uploaded file temporarily.
            with tempfile.NamedTemporaryFile(
                mode="wb",
                suffix=".fasta",
                delete=False
            ) as temp_file:
                for chunk in uploaded_file.chunks():
                    temp_file.write(chunk)

                temp_path = temp_file.name

            try:
                genome_data = parse_fasta(temp_path)

                if not genome_data:
                    raise ValueError("The FASTA file contains no sequences.")

                all_results = []
                total_bases = 0

                start_time = time.perf_counter()

                for chromosome, sequence in genome_data.items():
                    total_bases += len(sequence)

                    results = search_both_strands(
                        guide,
                        sequence,
                        max_mismatches=max_mismatches
                    )

                    for result in results:
                        result["chromosome"] = chromosome
                        all_results.append(result)

                execution_time = time.perf_counter() - start_time

            finally:
                os.remove(temp_path)

            context["analysis"] = {
                "total_records": len(genome_data),
                "total_bases": total_bases,
                "total_hits": len(all_results),
                "execution_time": round(execution_time, 6),
                "guide": guide,
                "max_mismatches": max_mismatches,
                "results": all_results,
                "status": "Completed",
                "mode": "CPU"
            }

        except (ValueError, OSError) as error:
            context["error"] = str(error)

    return render(request, "home.html", context)
