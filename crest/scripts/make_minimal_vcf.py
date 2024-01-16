# Create a minimal vcf for testing

from dataclasses import dataclass

header = '''##fileformat=VCFv4.2
##FILTER=<ID=LOW_TUMOR_VCN,Description="Germline variant has very low tumor variant copy number">
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allelic frequency calculated from read context counts as (Full + Partial + Core + Realigned + Alt) / Coverage">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=IMPACT,Number=10,Type=String,Description="Variant Impact [Gene, Transcript, CanonicalEffect, CanonicalCodingEffect, SpliceRegion, HgvsCodingImpact, HgvsProteinImpact, OtherReportableEffects, WorstCodingEffect, GenesAffected]">
##INFO=<ID=PURPLE_VCN,Number=1,Type=Float,Description="Purity adjusted variant copy number">
##INFO=<ID=DEVELOPER_COMMENT,Number=1,Type=String,Description="Developer Comment">
##contig=<ID=17,length=81195210>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	tumor_sample	ref_sample	rna_sample
'''


@dataclass
class Variant:
    chr: int
    pos: int
    ref: str
    alt: str
    filter: str  # e.g. "PASS", "LOW_TUMOR_VCN"
    gene: str
    ref_reads: int
    allele_reads: int
    total_reads: int
    comment: str

    def to_row(self) -> str:
        fields = (
            self.chr, self.pos, ".", self.ref, self.alt,
            500, # qual
            self.filter,
            self.info(),
            "GT:AD:AF:DP", # format
            self.tumor_sample(),
            self.ref_sample(),
            self.rna_sample(),
        )
        return '\t'.join((str(f) for f in fields)) + '\n'
 
    def info(self) -> str:
        if self.filter == "LOW_TUMOR_VCN":
            vcn = 0
        else:
            vcn = 1
        coding_effect = "NONE"
        worst_coding_effect = "NONE"
        comment = self.comment.replace(' ', '_')
        return f"IMPACT={self.gene},,,{coding_effect},,,,,{worst_coding_effect},1;PURPLE_VCN={vcn};DEVELOPER_COMMENT={comment}"

    def tumor_sample(self) -> str:
        return f"./.:0,100:1.0:100"

    def ref_sample(self) -> str:
        return f"1/1:0,30:1.0:30"

    def rna_sample(self) -> str:
        AD = f"{self.ref_reads},{self.allele_reads}"
        DP = f"{self.total_reads}"

        if self.total_reads  > 0:
            AF = self.allele_reads / self.total_reads  
        else:
            AF = 0.0
        return f"./.:{AD}:{AF}:{DP}"

def write_vcf(filename, data):
    with open(filename, "w") as f:
        f.write(header)

        for record in data:
            f.write(record.to_row())

if __name__ == '__main__':
    
    data = [
        Variant(17, 7579472, 'G', 'C', 'PASS', 'TP53', 48, 32, 80, "counted"),
        Variant(17, 7579473, 'G', 'C', 'LOW_TUMOR_VCN', 'TP53', 48, 0, 80, "not counted filter fail"),
        Variant(17, 7579474, 'G', 'C', 'PASS', '', 48, 32, 80, "not counted no gene impact"),
        Variant(17, 7579475, 'G', 'CC', 'PASS', 'TP53', 48, 32, 80, "not counted not a SNP"),
        Variant(17, 7579476, 'G', 'C', 'PASS', 'TP53', 48, 0, 80, "counted for total but not allele"),
        Variant(17, 7579477, 'G', 'C', 'PASS', 'TP53', 4, 1, 5, "not counted not enough total reads"),
    ]

    write_vcf("minimal.vcf", data)