
from geneweaver.parser.fasta_parser import (
    parse_fasta_in_chunks
)


def test_fasta_chunking(tmp_path):
    fasta = tmp_path / "sample.fasta"

    fasta.write_text(
        ">chromosome_1\nACGTACGTACGT\n"
    )

    chunks = list(
        parse_fasta_in_chunks(fasta, chunk_size=5)
    )

    assert len(chunks) == 3

    assert chunks[0]["sequence"] == "ACGTA"
    assert chunks[1]["sequence"] == "CGTAC"
    assert chunks[2]["sequence"] == "GT"

    assert chunks[0]["start"] == 0
    assert chunks[1]["start"] == 5
    assert chunks[2]["start"] == 10

    assert chunks[0]["header"] == "chromosome_1"
