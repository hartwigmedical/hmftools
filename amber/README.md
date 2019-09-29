# AMBER
AMBER is designed primarily to generate a tumor BAF file for use in PURPLE. 

AMBER locates heterozygous sites within the reference sample bam then calculates the allelic frequency of corresponding sites in the tumor bam. 
The Bioconductor copy number package is then used to generate pcf segments from the BAF file.

Additionally, AMBER is able to: 
  - detect evidence of contamination in the tumor from homozygous sites in the reference; and
  - facilitate sample matching by recording SNPs in the germline

## Installation

To install, download the latest compiled jar file from the [download links](#version-history-and-download-links). 
AMBER also requires a set of likely heterozygous sites (to determine BAF) and an optional set of SNP sites (to facilitate sample matching).  
Both are available to download from [HMFTools-Resources > Amber](https://resources.hartwigmedicalfoundation.nl/).

AMBER depends on the Bioconductor [copynumber](http://bioconductor.org/packages/release/bioc/html/copynumber.html) package for segmentation.
After installing [R](https://www.r-project.org/) or [RStudio](https://rstudio.com/), the copy number package can be added with the following R commands:
```
    library(BiocManager)
    install("copynumber")
```

AMBER requires Java 1.8+ and can be run with the minimum set of arguments as follows:

```
java -Xmx32G -cp amber.jar com.hartwig.hmftools.amber.AmberApplication \
   -reference COLO829R -reference_bam /run_dir/COLO829R.bam \ 
   -tumor COLO829T -tumor_bam /run_dir/COLO829T.bam \ 
   -output_dir /run_dir/amber/ \
   -threads 16 \
   -ref_genome /path/to/refGenome/refGenome.fasta \
   -bed /path/to/GermlineHetPon.hg19.bed.gz 
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
snp_bed| None | Locations to record SNPs in the germline
threads | 1 | Number of threads to use
min_mapping_quality | 1| Minimum mapping quality for an alignment to be used
min_base_quality | 13| Minimum quality for a base to be considered
min_depth_percent | 0.5 | Only include reference sites with read depth within min percentage of median reference read depth
max_depth_percent | 1.5 | Only include reference sites with read depth within max percentage of median reference read depth
min_het_af_percent | 0.4 | Minimum allelic frequency to be considered heterozygous
max_het_af_percent | 0.65 | Maximum allelic frequency to be considered heterozygous


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

## Output
File | Description
--- | ---
TUMOR.amber.baf.tsv | Tab separated values (TSV) containing reference and tumor BAF at each heterozygous site.
TUMOR.amber.baf.pcf | TSV of BAF segments using PCF algorithm.
TUMOR.amber.qc | Contains median tumor baf and QC status. FAIL may indicate contamination in sample. 
TUMOR.amber.baf.vcf.gz | Similar information as BAF file but in VCF format. 
TUMOR.amber.contamination.vcf.gz | Entry at each homozygous site in the reference and tumor.
REFERENCE.amber.snp.vcf.gz | Entry at each SNP location in the reference. 
 

# Version History and Download Links
- [2.5](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v2.5)
  - Fixed bug in contamination model if absolute zero contamination
- [2.4](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v2.4)
  - Added optional snp_bed parameter to output germline snps at specified locations
  - Changed file names and headers for better consistency with other HMF tools
- [2.3](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v2.3)
  - Gracefully handle contigs outside the ref genome. 
  - Fixed bug where TumorContamination file had two copies of tumor info rather than normal and tumor
  - CRAM support
- 2.2 
  - Fixed typo in header of TUMOR.amber.baf file.
- 2.1
  - Add statistical contamination check.
- 2.0
  - Read directly from BAMs without intermediary pileup step for significant performance improvements. 