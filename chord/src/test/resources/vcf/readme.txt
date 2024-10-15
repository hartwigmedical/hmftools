MINIMAL_SAMPLE.purple.somatic.vcf.gz: Variants were randomly subsampled using the commands:
zcat COLO829v003T.purple.somatic.vcf.gz | grep '#' | gzip -c > MINIMAL_SAMPLE.purple.somatic.vcf.gz &&
zcat COLO829v003T.purple.somatic.vcf.gz | grep -v '#' | grep PASS | shuf -n 2000 --random-source=COLO829v003T.purple.somatic.vcf.gz | sort -V | gzip -c >> MINIMAL_SAMPLE.purple.somatic.vcf.gz

MINIMAL_SAMPLE.purple.sv.vcf.gz: Variants were selected manually to cover DEL, DUP, INV, and BND

EMPTY_SAMPLE.purple.somatic.vcf.gz: Only VCF header from COLO829v003T.purple.somatic.vcf.gz

EMPTY_SAMPLE.purple.sv.vcf.gz: Only VCF header from COLO829v003T.purple.sv.vcf.gz