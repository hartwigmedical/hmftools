# AMBER (V2+)
AMBER is designed to generate a tumor BAF file for use in PURPLE. It locates heterozygous sites in the reference sample then calculates the allelic frequency of corresponding sites in the tumour.

## Inputs
AMBER requires:
 
1. Bams of both the normal and tumor sample. 
2. Reference genome file in fasta format
3. A bed file containing sites of likely heterozygous locations. The bed file used by HMF (GermlineHetPon.hg19.bed.gz) is available to download from [HMF-Pipeline-Resources.](https://resources.hartwigmedicalfoundation.nl) The sites were chosen by running the GATK HaplotypeCaller over 1700 germline samples and then selecting all SNP sites which are heterozygous in 800 to 900 of the samples. A HG38 bed file is also available.  

### Mandatory Arguments

Argument | Description 
---|---
-tumor | Name of the tumor sample
-tumor_bam | Path to tumor bam file
-reference | Name of the reference sample
-reference_bam | Path to reference bam file
-output_dir | Path to the output directory
-ref_genome | Path to the ref genome fasta file
-bed | Path to bed file containing likely heterozygous sites (refer above). Gz files supported.  

### Optional Arguments

Argument | Default | Description 
---|---|---
-threads | 1 | Number of threads to use
-min_mapping_quality | 1| Minimum mapping quality for an alignment to be used
-min_base_quality | 13| Minimum quality for a base to be considered
-min_depth_percent | 0.5 | Only include reference sites with read depth within min percentage of median reference read depth
-min_depth_percent | 1.5 | Only include reference sites with read depth within max percentage of median reference read depth
-min_het_af_percent | 0.4 | Minimum allelic frequency to be considered heterozygous
-max_het_af_percent | 0.65 | Maximum allelic frequency to be considered heterozygous

### Example Usage

```
java -cp amber.jar com.hartwig.hmftools.amber.AmberApplication \
   -tumor TUMOR \
   -tumor_bam /path/to/tumor/bam/TUMOR.bam \
   -reference REFERENCE \
   -reference_bam /path/to/reference/bam/REFERENCE.bam \
   -output_dir /run_dir/amber/ \
   -ref_genome /path/to/refGenome/refGenome.fasta \
   -bed /path/to/GermlineHetPon.hg19.bed.gz 
```

## Comparison to previous version
Calculating the BAF directly from the bams is functionalluy equivalent to the pileup method when using the following samtools arguments:

``` -A -B -x -Q 13 -q 1 -f /path/to/refGenome/refGenome.fasta ```


# AMBER (Pileup) - Deprecated

Note that this version is no longer supported. It is highly recommended you generate the BAFs direcly from the bams. 

## Prerequisites

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
java -cp amber.jar com.hartwig.hmftools.amber.pileup.AmberFromPileupApplication \ 
    -sample TUMOR \
    -output_dir /run_dir/amber \
    -reference /path/to/mpileup/REFERENCE.mpileup \
    -tumor /path/to/mpileup/TUMOR.mpileup
```

This will write output to `/run_dir/amber/TUMOR.amber.baf`

### Output

The output is a tab delimited file containing the Chromosome, Position and BAF of the tumor.