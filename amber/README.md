# AMBER

**A** **M**inipileup **B**AF  **E**stimato**r** is designed to generate a tumor BAF file for use in PURPLE.

Amber determines heterozygous locations of the reference sample then calculates the allelic frequency of corresponding locations in the tumour.

## Prerequsites

AMBER relies on mpileups of the reference and tumor samples sliced at likely heterozygous locations. Sambamba (or samtools) can be used to generate the mpileups.

Example generation:

```
export PATH=/path/to/samtools/:$PATH

sambamba mpileup -t 6 /
    -L /path/to/bed/HeterozygousLocations.bed /
    /path/to/bam/REFERENCE.bam /
    --samtools -q 1 -f /path/to/refGenome/refGenome.fasta /
    > /path/to/mpileup/REFERENCE.mpileup

sambamba mpileup -t 6 /
    -L /path/to/bed/HeterozygousLocations.bed /
    /path/to/bam/TUMOR.bam /
    --samtools -q 1 -f /path/to/refGenome/refGenome.fasta /
    > /path/to/mpileup/TUMOR.mpileup

```


## Usage

Argument | Default | Description
---|---|---
-sample | None | Name of tumor sample
-reference | None | Location of reference mpileup file
-tumor | None | Location of tumor mpileup file
-output_dir | None | Directory to write output
-min_depth_percent | 0.5 | Only include reference positions with read depth within min percentage of median
-min_depth_percent | 1.5 | Only include reference positions with read depth within max percentage of median
-min_het_af_percent | 0.4 | Minimum allelic frequency to be considered heterozygous
-max_het_af_percent | 0.65 | Maximum allelic frequency to be considered heterozygous

Arguments without default values are mandatory.

### Example Usage

```
java -jar amber.jar -sample TUMOR -output_dir /run_dir/amber -reference /path/to/mpileup/REFERENCE.mpileup -tumor /path/to/mpileup/TUMOR.mpileup
```

This will write output to `/run_dir/amber/TUMOR.amber.baf`

### Output

The output is a tab delimited file containing the Chromosome, Position and BAF of the tumor.