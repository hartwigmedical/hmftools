# COBALT
**Co**unt **ba**m **l**ines determines the read depth ratios of the supplied tumor and reference genomes. 

# Algorithm
The main steps are:
- reads are apportioned to windows to give raw depths and GC proportions per window
- filtering is applied to exclude windows that have unreliable data
- raw results are normalised to account for GC bias and other confounders
- a piecewise constant fitting algorithm is applied to segment the chromosomes to regions of constant copy number.

## Counting
### Window depth accumulation
Each read is assigned proportionally to the windows that it overlaps.
Reads of quality less than 10 (configurable) are ignored,
as are unmapped, duplicated, secondary and supplementary reads.
The raw count of bases from reads overlapping each window
is divided by the window length (1000) to give a raw depth.

### Window GC proportion
The proportion of Gs and Cs in bases of reads within a window is determined.
For low-depth windows (raw depth < 1.0), however, this value is overridden with the GC
value from the supplied GC profile data.

## Filtering
### Pseudo gene region filtering
Certain windows overlap pseudo genes on chromosomes 15 and 16 that produce artefactually high copy numbers.
These windows are excluded.

### Diploid region filtering
In tumor-only whole-genome mode, diploid regions can be specified, and windows that are outside of these regions are excluded.

### Target region filtering
In target region mode, only those windows that fall within the specified target regions are retained.

### GC Filtering
Windows with GC proportions outside the range [0.24,0.67] (configurable) are excluded.

### Mappability filtering
Windows with a mappability value less than 0.85 (configurable) are excluded.

## Normalisation
### GC normalisation
Each window is assigned to a bucket according to its GC proportion.
Median read depths  for allosome reads are calculated for each bucket.
The median read depths for the buckets are smoothed by averaging the value for each bucket with
those of its two neighbors. This produces a read depth normalisation factor for each bucket.
The read depth of each window is normalised by the factor for its associated GC bucket.
This produces a 'ratio' value which is the subject of all later normalisation steps
and which is the value used to detect copy number changes.

### Target region enrichment
In targeted mode, the ratio for each window is divided by a supplied enrichment factor.

### Median/mean normalisation
The median and mean read depths across all windows are calculated. The quotient of these
values is then applied as a normalisation factor across all windows.

### Diploid ratio normalisation
In whole-genome mode, windows are normalised to remove mega-base scale GC biases.
This normalization assumes that the median ratio of each 10Mb window (minimum 1Mb readable) 
should be diploid for autosomes and haploid for
male sex chromosomes, in addition to the following exceptions:

| Aberration                   | Chromosome       | Normalized Ratio   |
|------------------------------|------------------|--------------------|
| `MOSAIC_X`                   | X                | use median X ratio |
| `KLINEFELTER`                | X                | 1                  |
| `KLINEFELTER`                | Y                | 0.5                |
| `TRISOMY_[X,21,13,18,15,9P]` | X,21,13,18,15,9P | 1.5                |
| `TETRASOMY_9P`               | 9P               | 1.5                |


### Unity normalisation
The final normalisation step is to normalise so that in the output
file, the mean of the ratios for the windows is 1.0.

## Depth window consolidation
Sparse information in COBALT may cause a noisy fit for low pass WGS.
Therefore, we consolidate buckets to try to reach a median read depth of at least 8 per bucket. 
The ConsolidatedBucketSize is set to
```clamp(roundToOneSigDigit(80 / medianTumorReadCount, 10, 1000)```.
This formula allows consolidation into buckets of up to 1000 depth windows.
For standard WGS this should have no effect as medianTumorReadDepth >> 8.
We should never consolidate across regions of more than 3Mb (so never across centromere).
Consolidation is not used in targeted sequencing mode.

For the consolidated buckets, the mean GC ratio for both tumor and reference is calculated 
for the consolidated bucket and set to the centre window in the consolidated bucket.  
The other buckets are masked.

## Segmentation
Finally, the Bioconductor copy number package is used to generate segments from the ratio file.

All resource files for this tool and the WiGiTs pipeline are available for download via the [HMF Resource page](../pipeline/README_RESOURCES.md).

COBALT depends on the Bioconductor [copynumber](http://bioconductor.org/packages/release/bioc/html/copynumber.html) package for segmentation.
The R package [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html) is also used.
After installing [R](https://www.r-project.org/) or [RStudio](https://rstudio.com/), the required R packages can be added with the following R commands:
```
    library(BiocManager)
    install("copynumber")
    install("dplyr")
```

## Sequence and configuration of Cobalt stages
| Order | Step            | Type          | Parameters                     | Tumor, Targeted | Reference, Targeted | Tumor, Whole Genome | Reference, Whole Genome |
|-------|-----------------|---------------|--------------------------------|-----------------|---------------------|---------------------|-------------------------|
| 1     | Read depths     | counting      | `min_quality`                  | Y               | Y                   | Y                   | Y                       |
| 2     | GC proportion   | counting      | `gc_profile`                   | Y               | Y                   | Y                   | Y                       |
| 3     | Pseudo genes    | filtering     | built-in properties files      | Y               | Y                   | Y                   | Y                       |
| 3     | Mappability     | filtering     | `gc_profile`                   | Y               | Y                   | Y                   | Y                       |
| 3     | Diploid regions | filtering     | `tumor_only_diploid_bed`       | N               | N                   | Y                   | Y                       |
| 4     | Target regions  | filtering     | `target_region_norm_file`      | Y               | Y                   | N                   | N                       |
| 5     | GC              | filtering     | `gc_ratio_min`, `gc_ratio_max` | Y               | Y                   | Y                   | Y                       |
| 5     | GC              | normalisation | `gc_profile`                   | Y               | Y                   | Y                   | Y                       |
| 6     | Enrichment      | normalisation | `target_region_norm_file`      | Y               | Y                   | N                   | N                       |
| 7     | Median/mean     | normalisation |                                | N               | N                   | Y                   | Y                       |
| 8     | Low coverage    | consolidation |                                | N               | N                   | Y                   | Y                       |
| 9     | Diploid ratios  | normalisation |                                | N               | N                   | Y                   | Y                       |
| 10    | Unity           | normalisation |                                | Y               | Y                   | N                   | N                       |
| 11    | Copy number     | segmentation  | `pcf_gamma`, `skip_pcf_calc`   | Y               | Y                   | Y                   | Y                       |


# Running Cobalt
## Tumor-only Mode
In the absence of a reference bam and reference COBALT will  be run in tumor-only mode.

Without a means to determine which regions of the normal are diploid, a bed file specifying these
locations must be included with the `tumor-only-diploid-bed` parameter.
A 37 bed file (DiploidRegions.37.bed.gz) and 38 equivalent are available to download
from [HMF-Pipeline-Resources](https://resources.hartwigmedicalfoundation.nl). To create this bed file we examined the COBALT output of 100 samples.  
We considered each 1000 base region to be diploid if 50% or more of the samples
were diploid (0.85 >= referenceGCDiploidRatio <= 1.15 ) at this point.

As no reference data is supplied, COBALT does not try to determine gender or any chromosomal aberrations.
No reference pcf file will be created. The output reference ratios will be 1 or -1 on all chromosomes even
if they are allosomes.  Downstream, PURPLE will adjust the allosome ratios according to the AMBER gender.

## Germline Only Mode
In the absence of a tumor bam and tumor COBALT will  be run in germline mode.
Counts and ratios are only calculated and fitted for the reference sample.

## Targeted Mode
COBALT may be run on targeted data. For more information on how to run hmftools in targeted mode please see [here](https://github.com/hartwigmedical/hmftools/blob/master/pipeline/README_TARGETED.md)

## Arguments
COBALT requires Java 17+ and for tumor-only mode can be run with the minimum set of arguments as follows:
```
java -jar -Xmx8G cobalt.jar \
    -tumor SAMPLE_ID \
    -tumor_bam /sample_data/SAMPLE_ID.bam \ 
    -output_dir /sample_data/ \ 
    -threads 10 \ 
    -ref_genome_version 37
    -gc_profile /ref_data/GC_profile.1000bp.37.cnp
```

For in tumor-reference mode the reference sample must be supplied:
```
java -jar -Xmx8G cobalt.jar \
    -reference SAMPLE_ID_R \
    -reference_bam /sample_data/SAMPLE_ID_R.bam \ 
    -tumor SAMPLE_ID \
    -tumor_bam /sample_data/SAMPLE_ID.bam \ 
    -output_dir /sample_data/ \ 
    -threads 10 \ 
    -ref_genome_version 37
    -gc_profile /ref_data/GC_profile.1000bp.37.cnp
```

### Mandatory Arguments
| Argument           | Description                                                                               |
|--------------------|-------------------------------------------------------------------------------------------|
| reference          | Name of the reference sample (not needed if tumor supplied)                               |
| reference_bam      | Path to reference BAM file                                                                |
| tumor              | Name of tumor sample (not needed if reference supplied)                                   |
| tumor_bam          | Path to tumor BAM file                                                                    |
| output_dir         | Path to the output directory. This directory will be created if it does not already exist |
| gc_profile         | Path to GC profile                                                                        |
| ref_genome_version | One of `37` or `38`                                                                       |

A compressed copy of the GC Profile file used by HMF (GC_profile.1000bp.37.cnp) is available to 
download from [HMF-Pipeline-Resources](https://resources.hartwigmedicalfoundation.nl). This file contains 5 columns for each 1kb window 
of the genome {chromosome,position,GC Proportion,Non N Proportion,Average Mappability}.  
A 38 equivalent is also available. Please note the downloaded file must be un-compressed before use. 

COBALT supports both BAM and CRAM file formats. If using CRAM, the ref_genome argument must be included.

### Optional Arguments
| Argument               | Default | Description                                                 |
|------------------------|---------|-------------------------------------------------------------|
| threads                | 4       | Number of threads to use                                    |
| min_quality            | 10      | Min quality                                                 |
| ref_genome             | None    | Path to the reference genome fasta file if using CRAM files |
| validation_stringency  | STRICT  | SAM validation strategy: STRICT, SILENT, LENIENT            |
| tumor_only_diploid_bed | NA      | Bed file of diploid regions of the genome                   |
| pcf_gamma              | 100     | Gamma value for use in R copy_number pcf function           |
| target_region          | None    | Target region TSV file for use in targeted mode.            |

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

# Output
The following tab delimited files are written:

`/run_dir/cobalt/TUMOR.cobalt.ratio.tsv.gz`

`/run_dir/cobalt/TUMOR.cobalt.ratio.pcf`

`/run_dir/cobalt/REFERENCE.cobalt.ratio.pcf`

TUMOR.cobalt.ratio.tsv.gz contains the counts and ratios of the reference and tumor
and the GC proportions:

| chromosome | position | referenceReadDepth | tumorReadDepth | referenceGCRatio | tumorGCRatio | referenceGCDiploidRatio | referenceGCContent | tumorGCContent |
|------------|----------|--------------------|----------------|------------------|--------------|-------------------------|--------------------|----------------|
| 1          | 4000001  | 37.367             | 143.34         | 0.8827           | 0.9754       | 0.8904                  | 0.6                | 0.6022         |
| 1          | 4001001  | 41.327             | 119.101        | 1.0626           | 0.8792       | 1.0719                  | 0.5433             | 0.5391         |
| 1          | 4002001  | 35.33              | 127.256        | 0.8346           | 0.8802       | 0.8418                  | 0.5994             | 0.5941         |
| 1          | 4003001  | 35.928             | 111.033        | 0.9737           | 0.8496       | 0.9822                  | 0.4767             | 0.4753         |
| 1          | 4004001  | 36.938             | 122.323        | 0.8725           | 0.8324       | 0.8801                  | 0.5973             | 0.5985         |

TUMOR.cobalt.ratio.pcf and REFERENCE.cobalt.ratio.pcf contain the segmented regions determined from the ratios.


# Version History and Download Links
- [1.14.1](https://github.com/hartwigmedical/hmftools/releases/tag/cobalt-v1.14.1)
  - Fix crash bug in the bucket consolidation.
- [1.14](https://github.com/hartwigmedical/hmftools/releases/tag/cobalt-v1.14)
  - Automatically consolidating buckets if mean coverage <= 50.
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
