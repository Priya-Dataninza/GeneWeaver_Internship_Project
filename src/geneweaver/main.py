
from geneweaver.parser.fasta_parser import parse_fasta
from geneweaver.cpu.baseline import search_off_targets


def main():
    guide = "ACGTACGT"

    genome_data = parse_fasta("data/sample_genome.fasta")

    for chromosome, sequence in genome_data.items():
        results = search_off_targets(
            guide,
            sequence,
            max_mismatches=2
        )

        print(f"\n{chromosome}")
        print(f"Candidates found: {len(results)}")

        for result in results:
            print(result)


if __name__ == "__main__":
    main()
