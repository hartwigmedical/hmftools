# COBALT

**Co**unt **ba**m **l**ines determines the read depth ratios of the supplied tumor and reference genomes. 

COBALT starts with the raw read counts per 1,000 base window for both normal and tumor samples by counting the number of alignment starts 
in the respective bam files with a mapping quality score of at least 10 that is neither unmapped, duplicated, secondary, nor supplementary. 
Windows with a GC content less than 0.2 or greater than 0.6 or with an average mappability below 0.85 are excluded from further analysis.

Next we apply a GC normalization to calculate the read ratios. 
To do this we divide the read count of each window by the median read count of all windows sharing the same GC content then normalise further to the 
ratio of the median to mean read count of all windows. 

Post GC normalization, COBALT is able to detect the following germline chromosomal aberrations from the reference ratio:

Aberration | Gender | Ratio Criteria
---|---|---
`MOSAIC_X` | FEMALE| X ratio < min(0.8, minAutosomeMedianDepthRatio*)
`KLINEFELTER` (XXY) | MALE | X ratio >= 0.65
`TRISOMY_[X,21,13,18,15]` | BOTH | chromosome ratio >= 1.4

*By checking against autosomes we rule out very high GC bias in the reference.  

The reference sample ratios have a further ‘diploid’ normalization applied to them to remove megabase scale GC biases. 
This normalization assumes that the median ratio of each 10Mb window (minimum 1Mb readable) should be diploid for autosomes and haploid for 
male sex chromosomes in addition to the following exceptions:

Aberration | Chromosome | Normalized Ratio
---|---|---
`MOSAIC_X` | X | use median X ratio
`KLINEFELTER` | X | 1
`KLINEFELTER` | Y | 0.5
`TRISOMY_[X,21,13,18,15]` | X,21,13,18,15 | 1.5

Finally, the Bioconductor copy number package is used to generate segments from the ratio file.

## Installation

To install, download the latest compiled jar file from the [download links](#version-history-and-download-links) and the appropriate GC profile from [HMFTools-Resources > Cobalt](https://resources.hartwigmedicalfoundation.nl/).

COBALT depends on the Bioconductor [copynumber](http://bioconductor.org/packages/release/bioc/html/copynumber.html) package for segmentation.
After installing [R](https://www.r-project.org/) or [RStudio](https://rstudio.com/), the copy number package can be added with the following R commands:
```
    library(BiocManager)
    install("copynumber")
```

COBALT requires Java 1.8+ and can be run with the minimum set of arguments as follows:

```
java -cp -Xmx8G cobalt.jar com.hartwig.hmftools.cobalt.CobaltApplication \
    -reference COLO829R -reference_bam /run_dir/COLO829R.bam \ 
    -tumor COLO829T -tumor_bam /run_dir/COLO829T.bam \ 
    -output_dir /run_dir/cobalt \ 
    -threads 16 \ 
    -gc_profile /path/to/GC_profile.1000bp.37.cnp
```

## Mandatory Arguments

Argument  | Description
---|---
reference | Name of the reference sample
reference_bam | Path to reference BAM file
tumor | Name of tumor sample
tumor_bam | Path to tumor BAM file
output_dir | Path to the output directory. This directory will be created if it does not already exist
gc_profile | Path to GC profile 

A compressed copy of the GC Profile file used by HMF (GC_profile.1000bp.37.cnp) is available to download from [HMF-Pipeline-Resources](https://resources.hartwigmedicalfoundation.nl). 
A 38 equivalent is also available. Please note the downloaded file must be un-compressed before use. 

COBALT supports both BAM and CRAM file formats. If using CRAM, the ref_genome argument must be included.

## Optional Arguments

Argument | Default | Description 
---|---|---
threads | 4 | Number of threads to use
min_quality | 10 | Min quality
ref_genome | None | Path to the reference genome fasta file if using CRAM files
validation_stringency | STRICT | SAM validation strategy: STRICT, SILENT, LENIENT
tumor_only | NA | Set to tumor only mode
tumor_only_diploid_bed | NA | Bed file of diploid regions of the genome

## Tumor Only Mode
In the absence of a reference bam, COBALT can be put into tumor only mode with the `tumor_only` flag. 
In this mode the `reference` and `reference_bam` parameters are no longer valid.

As no reference data is supplied, COBALT does not try to determine gender or any chromosomal aberrations. 
The output reference ratios will be 1 or -1 on all chromosomes even if they are allosomes. 
Downstream, PURPLE will adjust the allosome ratios according to the AMBER gender. 
A file called `DIPLOID.cobalt.ratio.pcf' will be created in lieu of a reference PCF file.

Without a means to determine which regions of the normal are diploid, a bed file specifying these locations must be included with the `tumor-only-diploid-bed` parameter. 
A 37 bed file (DiploidRegions.37.bed.gz) and 38 equivalent are available to download from [HMF-Pipeline-Resources](https://resources.hartwigmedicalfoundation.nl).

To create this bed file we examined the COBALT output of 100 samples. 
We considered each 1000 base region to be diploid if 50% or more of the samples were diploid (0.85 >= referenceGCDiploidRatio <= 1.15 ) at this point. 

## Performance Characteristics
Performance numbers were taken from a 72 core machine using COLO829 data with an average read depth of 35 and 93 in the normal and tumor respectively. 
Elapsed time is measured in minutes. 
CPU time is minutes spent in user mode. 
Peak memory is measure in gigabytes.


Threads | Elapsed Time| CPU Time | Peak Mem
---|---|---|---
1 | 111 | 122 | 3.85
8 | 17 | 127 | 4.49
16 | 10 | 139 | 4.58 
32 | 11 | 184 | 4.33
48 | 10 | 153 | 4.35


## Output
The following tab delimited files are written:

`/run_dir/cobalt/TUMOR.chr.len`

`/run_dir/cobalt/TUMOR.cobalt.ratio.tsv`

`/run_dir/cobalt/TUMOR.cobalt.ratio.pcf`

`/run_dir/cobalt/REFERENCE.cobalt.ratio.pcf`

TUMOR.cobalt.ratio.tsv contains the counts and ratios of the reference and tumor:

Chromosome | Position | ReferenceReadCount | TumorReadCount | ReferenceGCRatio | TumorGCRatio | ReferenceGCDiploidRatio
---|---|---|---|---|---|---
1|4000001|204|504|0.8803|0.855|0.8982
1|4001001|203|570|0.8429|0.9149|0.86
1|4002001|155|473|0.6463|0.7654|0.6594
1|4003001|260|566|1.098|0.9328|1.1203
1|4004001|256|550|1.1144|0.9428|1.1371

TUMOR.cobalt.ratio.pcf and REFERENCE.cobalt.ratio.pcf contain the segmented regions determined from the ratios.

## Migration to 1.9

As germline aberrations only effect the final normalization it is possible to migrate existing COBALT output from versions 1.4 to 1.8 to 1.9 without having to re-examine the bams using the following command:

```
java -Xmx8G -Xms4G -cp ${cobalt_jar} com.hartwig.hmftools.cobalt.CountBamLinesMigration \
    -reference COLO829R -tumor COLO829T \
    -gc_profile /path/to/GC_profile.1000bp.37.cnp \
    -input_dir /path/to/existing/cobalt/data/ \
    -output_dir /path/to/new/cobalt/data/ 
```

## Version History and Download Links
- [1.11](https://github.com/hartwigmedical/hmftools/releases/tag/cobalt-v1.11)
  - Tumor only mode
- [1.10](https://github.com/hartwigmedical/hmftools/releases/tag/cobalt-v1.10)
  - Re-added support for cancel panel integration test 
- [1.9](https://github.com/hartwigmedical/hmftools/releases/tag/cobalt-v1.9)
  - Alert user that gc_profile should be un-compressed before use
  - Add support for XXY, XYY, Female Mosaic X, and Trisomy 13,15,18,21,X
- [1.8](https://github.com/hartwigmedical/hmftools/releases/tag/cobalt-v1.8)
  - Added `validation_stringency` parameter.
  - Added explicit `stringsAsFactors = T` to R script
- [1.7](https://github.com/hartwigmedical/hmftools/releases/tag/cobalt-v1.7)
  - Exit gracefully on exceptions
  - Changed file names and headers for better consistency with other HMF tools
- [1.6](https://github.com/hartwigmedical/hmftools/releases/tag/cobalt-v1.6)
  - CRAM support
- 1.5
  - Support for ref genome 38
- 1.4
  - Support for Klinefelter syndrome