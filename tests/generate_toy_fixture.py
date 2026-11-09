"""Rebuild the tiny synthetic BAM and indexed VCFs in toy_data/."""

from pathlib import Path
import tempfile

import pysam


POSITIONS = (100, 120, 140, 160, 180)
# Each character specifies REF (0) or ALT (1) at the five somatic sites.
PATTERNS = [("00000", 25), ("11000", 5), ("11100", 5), ("11110", 10), ("11001", 5)]


def main():
    output = Path(__file__).resolve().parents[1] / "toy_data"
    output.mkdir(exist_ok=True)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    read_number = 0
    with pysam.AlignmentFile(str(output / "reads.bam"), "wb", header=header) as bam:
        for pattern, count in PATTERNS:
            sequence = list("A" * 250)
            for position, allele in zip(POSITIONS, pattern):
                if allele == "1":
                    sequence[position - 1] = "G"
            if "1" in pattern:
                for position in (20, 40):
                    sequence[position - 1] = "G"
            for _ in range(count):
                read_number += 1
                read = pysam.AlignedSegment()
                read.query_name = f"toy{read_number:02d}"
                read.query_sequence = "".join(sequence)
                read.query_qualities = pysam.qualitystring_to_array("I" * 250)
                read.reference_id = 0
                read.reference_start = 0
                read.mapping_quality = 60
                read.cigarstring = "250M"
                bam.write(read)
    pysam.index(str(output / "reads.bam"))

    with tempfile.TemporaryDirectory() as temporary:
        for name, positions, germline in [("somatic", POSITIONS, False), ("germline", (20, 40), True)]:
            text = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
            if germline:
                text += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            text += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
            text += "\tFORMAT\tNORMAL\n" if germline else "\n"
            for position in positions:
                text += f"chr1\t{position}\t.\tA\tG\t.\tPASS\t."
                text += "\tGT\t1|0\n" if germline else "\n"
            plain = Path(temporary) / f"{name}.vcf"
            plain.write_text(text)
            compressed = str(output / f"{name}.vcf.gz")
            pysam.tabix_compress(str(plain), compressed, force=True)
            pysam.tabix_index(compressed, preset="vcf", force=True)
    print(f"Wrote {read_number} reads and indexed VCFs to {output}")


if __name__ == "__main__":
    main()
