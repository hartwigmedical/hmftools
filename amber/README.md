# AMBER
AMBER is designed to generate a tumor BAF file for use in PURPLE. 

Looking directly at the BAM files, AMBER locates heterozygous sites within the reference sample then calculates the allelic frequency of corresponding sites in the tumor. 

AMBER also locates homozygous sites in the reference to detect evidence of contamination in the tumor.

Finally, the Bioconductor copy number package is used to generate pcf segments from the BAF file.

Prior versions of AMBER relied on mpileups of the tumor and reference rather than the BAMs themselves. 
This method is deprecated but still available in the jar file as AmberFromPileupApplication. 
See below for more details.  

## R Dependencies
Segmentation is done with the Bioconductor [copynumber](http://bioconductor.org/packages/release/bioc/html/copynumber.html) package.

This can be installed in R (3.5+) with the following commands:
```
    library(BiocManager)
    install("copynumber")
```

## Mandatory Arguments

Argument | Description 
---|---
reference | Name of the reference sample
reference_bam | Path to indexed reference BAM file
tumor | Name of the tumor sample
tumor_bam | Path to indexed tumor BAM file
output_dir | Path to the output directory. This directory will be created if it does not already exist.
ref_genome | Path to the reference genome fasta file
bed | Path to bed file containing likely heterozygous sites (see below). Gz files supported.  

The bed file used by HMF (GermlineHetPon.hg19.bed.gz) is available to download from [HMF-Pipeline-Resources](https://resources.hartwigmedicalfoundation.nl). 
The sites were chosen by running the GATK HaplotypeCaller over 1700 germline samples and then selecting all SNP sites which are heterozygous in 800 to 900 of the samples. 
The 1.3 million sites provided in this file typically result in 450k+ BAF points. A HG38 equivalent is also available.

AMBER supports both BAM and CRAM file formats. 

## Optional Arguments

Argument | Default | Description 
---|---|---
threads | 1 | Number of threads to use
min_mapping_quality | 1| Minimum mapping quality for an alignment to be used
min_base_quality | 13| Minimum quality for a base to be considered
min_depth_percent | 0.5 | Only include reference sites with read depth within min percentage of median reference read depth
max_depth_percent | 1.5 | Only include reference sites with read depth within max percentage of median reference read depth
min_het_af_percent | 0.4 | Minimum allelic frequency to be considered heterozygous
max_het_af_percent | 0.65 | Maximum allelic frequency to be considered heterozygous

## Example Usage

```
java -Xmx32G -cp amber.jar com.hartwig.hmftools.amber.AmberApplication \
   -reference COLO829R -reference_bam /run_dir/COLO829R.bam \ 
   -tumor COLO829T -tumor_bam /run_dir/COLO829T.bam \ 
   -output_dir /run_dir/amber/ \
   -threads 16 \
   -ref_genome /path/to/refGenome/refGenome.fasta \
   -bed /path/to/GermlineHetPon.hg19.bed.gz 
```


## Performance Characteristics
Performance numbers were taken from a 72 core machine using COLO829 data with an average read depth of 35 and 93 in the normal and tumor respectively. 
Elapsed time is measured in minutes. 
CPU time is minutes spent in user mode. 
Peak memory is measure in gigabytes.

Reading directly from the bam, AMBER has the following characteristics:

Amber Method | Threads | Elapsed Time| CPU Time | Peak Mem
---|---|---|---|---
Bam | 1 | 144 | 230 | 15.04
Bam | 8 | 22 | 164 | 18.40
Bam | 16 | 12 | 164 | 21.00
Bam | 32 | 8 | 170 | 21.60
Bam | 48 | 7 | 199 | 21.43
Bam | 64 | 6 | 221 | 21.78

For comparison, the deprecated pileup method has the following characteristics:

Amber Method | Threads | Elapsed Time| CPU Time 
---|---|---|---
Pileup | 1 | 180 | 180 
Pileup | 64 | 44 | 272 

## Output
File | Description
--- | ---
TUMOR.amber.baf | Tab separated values (TSV) containing reference and tumor BAF at each heterozygous site.
TUMOR.amber.baf.pcf | TSV of BAF segments using PCF algorithm.
TUMOR.amber.qc | Contains median tumor baf and QC status. FAIL may indicate contamination in sample. 
TUMOR.amber.vcf.gz | Similar information as BAF file but in VCF format. This file is not used by PURPLE.
 
## Comparison to pileup method
Calculating the BAF directly from the bams is functionally equivalent to the pileup method when using the following samtools arguments:

``` -A -B -x -Q 13 -q 1 -f /path/to/refGenome/refGenome.fasta ```

# AMBER From Pileup - Deprecated

This version of AMBER relies on mpileups of the reference and tumor samples sliced at likely heterozygous locations. Sambamba (or samtools) can be used to generate the mpileups. 

Example generation:

```
export PATH=/path/to/samtools/:$PATH

sambamba mpileup -t 6 /
    -L /path/to/bed/HeterozygousLocations.bed /
    /path/to/bam/REFERENCE.bam /
    --samtools -A -B -x -Q 13 -q 1 -f /path/to/refGenome/refGenome.fasta /
    > /path/to/mpileup/REFERENCE.mpileup

sambamba mpileup -t 6 /
    -L /path/to/bed/HeterozygousLocations.bed /
    /path/to/bam/TUMOR.bam /
    --samtools -A -B -x -Q 13 -q 1 -f /path/to/refGenome/refGenome.fasta /
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

## Example Usage

Please note the application name has changed from earlier versions.

```
java -cp amber.jar com.hartwig.hmftools.amber.pileup.AmberFromPileupApplication \ 
    -sample TUMOR \
    -output_dir /run_dir/amber \
    -reference /path/to/mpileup/REFERENCE.mpileup \
    -tumor /path/to/mpileup/TUMOR.mpileup
```

## Version History
- [2.3](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v2-3)
  - Gracefully handle contigs outside the ref genome. 
  - Fixed bug where TumorContamination file had two copies of tumor info rather than normal and tumor
  - CRAM support
- 2.2 
  - Fixed typo in header of TUMOR.amber.baf file.
- 2.1
  - Add statistical contamination check.
- 2.0
  - Read directly from BAMs without intermidiary pileup step for significant performance improvements. 