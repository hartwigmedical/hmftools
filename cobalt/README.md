# COBALT

**Co**unt **ba**m **l**ines determines the read depth ratios of the supplied tumor and reference genomes. 

COBALT starts with the raw read counts per 1,000 base window for both normal and tumor samples by counting the number of alignment starts 
in the respective bam files with a mapping quality score of at least 10 that is neither unmapped, duplicated, secondary, nor supplementary. 
Windows with a GC content less than 0.2 or greater than 0.6 or with an average mappability below 0.85 are excluded from further analysis.

Next we apply a GC normalization to calculate the read ratios. To do this we divide the read count of each window by the median read count of all windows sharing the same GC content then normalise further to the ratio of the median to mean read count of all windows.    For some targeted region analyses the median read count may be very low due to low numbers of off target reads, and the use of discrete integers may restrict resolution.  We therefore modify the median with the following formula to improve the estimate:

```
modifiedMedian = median - 0.5 + (0.5 - proportion of depth windows with read count < median) / (proportion of depth windows  with read count = median)
```

Post GC normalization, COBALT is able to detect the following germline chromosomal aberrations from the reference ratio:

| Aberration                | Gender | Ratio Criteria                                   |
|---------------------------|--------|--------------------------------------------------|
| `MOSAIC_X`                | FEMALE | X ratio < min(0.8, minAutosomeMedianDepthRatio*) |
| `KLINEFELTER` (XXY)       | MALE   | X ratio >= 0.65                                  |
| `TRISOMY_[X,21,13,18,15]` | BOTH   | chromosome ratio >= 1.4                          |

*By checking against autosomes we rule out very high GC bias in the reference.  

The reference sample ratios have a further ‘diploid’ normalization applied to them to remove megabase scale GC biases. 
This normalization assumes that the median ratio of each 10Mb window (minimum 1Mb readable) should be diploid for autosomes and haploid for 
male sex chromosomes in addition to the following exceptions:

| Aberration                | Chromosome    | Normalized Ratio   |
|---------------------------|---------------|--------------------|
| `MOSAIC_X`                | X             | use median X ratio |
| `KLINEFELTER`             | X             | 1                  |
| `KLINEFELTER`             | Y             | 0.5                |
| `TRISOMY_[X,21,13,18,15]` | X,21,13,18,15 | 1.5                |

Finally, the Bioconductor copy number package is used to generate segments from the ratio file.

## Installation

To install, download the latest compiled jar file from the [download links](#version-history-and-download-links) and the appropriate GC profile from [HMFTools-Resources > DNA Pipeline](https://storage.googleapis.com/hmf-public/HMFtools-Resources/dna_pipeline/).

COBALT depends on the Bioconductor [copynumber](http://bioconductor.org/packages/release/bioc/html/copynumber.html) package for segmentation.
After installing [R](https://www.r-project.org/) or [RStudio](https://rstudio.com/), the copy number package can be added with the following R commands:
```
    library(BiocManager)
    install("copynumber")
```

COBALT requires Java 11+ and can be run with the minimum set of arguments as follows:

```
java -jar -Xmx8G cobalt.jar \
    -reference COLO829R -reference_bam /run_dir/COLO829R.bam \ 
    -tumor COLO829T -tumor_bam /run_dir/COLO829T.bam \ 
    -output_dir /run_dir/cobalt \ 
    -threads 16 \ 
    -gc_profile /path/to/GC_profile.1000bp.37.cnp
```

## Mandatory Arguments

| Argument      | Description                                                                               |
|---------------|-------------------------------------------------------------------------------------------|
| reference     | Name of the reference sample                                                              |
| reference_bam | Path to reference BAM file                                                                |
| tumor         | Name of tumor sample                                                                      |
| tumor_bam     | Path to tumor BAM file                                                                    |
| output_dir    | Path to the output directory. This directory will be created if it does not already exist |
| gc_profile    | Path to GC profile                                                                        |

A compressed copy of the GC Profile file used by HMF (GC_profile.1000bp.37.cnp) is available to download from [HMF-Pipeline-Resources](https://resources.hartwigmedicalfoundation.nl). 
A 38 equivalent is also available. Please note the downloaded file must be un-compressed before use. 

COBALT supports both BAM and CRAM file formats. If using CRAM, the ref_genome argument must be included.

## Optional Arguments

| Argument               | Default | Description                                                 |
|------------------------|---------|-------------------------------------------------------------|
| threads                | 4       | Number of threads to use                                    |
| min_quality            | 10      | Min quality                                                 |
| ref_genome             | None    | Path to the reference genome fasta file if using CRAM files |
| validation_stringency  | STRICT  | SAM validation strategy: STRICT, SILENT, LENIENT            |
| tumor_only_diploid_bed | NA      | Bed file of diploid regions of the genome                   |
| pcf_gamma              | 100     | Gamma value for use in R copy_number pcf function           |
| target_region          | None    | Target region TSV file for use in targeted mode.            |

## Tumor Only Mode
In the absence of a reference bam and reference COBALT will  be run in tumor_only mode.    

Without a means to determine which regions of the normal are diploid, a bed file specifying these locations must be included with the `tumor-only-diploid-bed` parameter. A 37 bed file (DiploidRegions.37.bed.gz) and 38 equivalent are available to download from [HMF-Pipeline-Resources](https://resources.hartwigmedicalfoundation.nl). To create this bed file we examined the COBALT output of 100 samples.  We considered each 1000 base region to be diploid if 50% or more of the samples were diploid (0.85 >= referenceGCDiploidRatio <= 1.15 ) at this point. 

As no reference data is supplied, COBALT does not try to determine gender or any chromosomal aberrations.  No reference pcf file will be created. The output reference ratios will be 1 or -1 on all chromosomes even if they are allosomes.  Downstream, PURPLE will adjust the allosome ratios according to the AMBER gender. 

## Germline Only Mode

In the absence of a tumor bam and tumor COBALT will  be run in germline mode.   Counts and ratios are only calculated and fitted for the reference sample.

## Targeted Mode

COBALT may be run on targeted data.   For more information on how to run hmftools in targeted mode please see [here](https://github.com/hartwigmedical/hmftools/blob/master/README_TARGETED.md)

## Performance Characteristics
Performance numbers were taken from a 72 core machine using COLO829 data with an average read depth of 35 and 93 in the normal and tumor respectively. 
Elapsed time is measured in minutes. 
CPU time is minutes spent in user mode. 
Peak memory is measure in gigabytes.


| Threads | Elapsed Time | CPU Time | Peak Mem |
|---------|--------------|----------|----------|
| 1       | 111          | 122      | 3.85     |
| 8       | 17           | 127      | 4.49     |
| 16      | 10           | 139      | 4.58     |
| 32      | 11           | 184      | 4.33     |
| 48      | 10           | 153      | 4.35     |

## Output
The following tab delimited files are written:

`/run_dir/cobalt/TUMOR.cobalt.ratio.tsv.gz`

`/run_dir/cobalt/TUMOR.cobalt.ratio.pcf`

`/run_dir/cobalt/REFERENCE.cobalt.ratio.pcf`

TUMOR.cobalt.ratio.tsv.gz contains the counts and ratios of the reference and tumor:

| Chromosome | Position | ReferenceReadCount | TumorReadCount | ReferenceGCRatio | TumorGCRatio | ReferenceGCDiploidRatio |
|------------|----------|--------------------|----------------|------------------|--------------|-------------------------|
| 1          | 4000001  | 204                | 504            | 0.8803           | 0.855        | 0.8982                  |
| 1          | 4001001  | 203                | 570            | 0.8429           | 0.9149       | 0.86                    |
| 1          | 4002001  | 155                | 473            | 0.6463           | 0.7654       | 0.6594                  |
| 1          | 4003001  | 260                | 566            | 1.098            | 0.9328       | 1.1203                  |
| 1          | 4004001  | 256                | 550            | 1.1144           | 0.9428       | 1.1371                  |

TUMOR.cobalt.ratio.pcf and REFERENCE.cobalt.ratio.pcf contain the segmented regions determined from the ratios.


## Version History and Download Links
- [1.13](https://github.com/hartwigmedical/hmftools/releases/tag/cobalt-v1.13)
  - Added support for germline only mode.
  - Added support for targeted mode. Activated when run with `-target_region` argument.
  - Keeps 0 read counts for regions that have no read, instead of dropping those regions from output.
  - Added new argument `-pcf_gamma` for overriding PCF gamma value.
  - Remove some redundant output files.
- [1.12](https://github.com/hartwigmedical/hmftools/releases/tag/cobalt-v1.12)
  - Added workaround for R copy_number module pcf function bug 
- [1.11](https://github.com/hartwigmedical/hmftools/releases/tag/cobalt-v1.11)
  - Tumor only mode
- [1.10](https://github.com/hartwigmedical/hmftools/releases/tag/cobalt-v1.10)
  - Re-added support for cancel panel integration test 
- [1.9](https://github.com/hartwigmedical/hmftools/releases/tag/cobalt-v1.9)
  - Alert user that gc_profile should be un-compressed before use
  - Add support for XXY, XYY, Female Mosaic X, and Trisomy 13,15,18,21,X
