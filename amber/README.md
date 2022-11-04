# AMBER
AMBER is designed to generate a tumor BAF file for use in PURPLE from a provided VCF of likely heterozygous SNP sites.

When using paired reference/tumor bams, AMBER confirms these sites as heterozygous in the reference sample bam then calculates the allelic frequency of corresponding sites in the tumor bam. 
In tumor only mode, all provided sites are examined in the tumor with additional filtering then applied. 
 
The Bioconductor copy number package is then used to generate pcf segments from the BAF file.

When using paired reference/tumor data, AMBER is also able to: 
  - detect evidence of contamination in the tumor from homozygous sites in the reference; and
  - facilitate sample matching by recording SNPs in the germline
  - identify long regions of homozygosty and consanguinity

## Installation

To install, download the latest compiled jar file from the [download links](#version-history-and-download-links). 

Ref genome versions 37 and 38 of the likely heterozygous sites are available to download from [HMFTools-Resources > DNA Pipeline > copy_number](https://storage.googleapis.com/hmf-public/HMFtools-Resources/dna_pipeline/).

The Bioconductor [copynumber](http://bioconductor.org/packages/release/bioc/html/copynumber.html) package is required for segmentation.
After installing [R](https://www.r-project.org/) or [RStudio](https://rstudio.com/), the copy number package can be added with the following R commands:
```
    library(BiocManager)
    install("copynumber")
```

AMBER requires Java 11+ to be installed.

## Paired Normal/Tumor Mode
This is the default and recommended mode.

### Mandatory Arguments

| Argument      | Description                                                                                |
|---------------|--------------------------------------------------------------------------------------------|
| reference     | Name of the reference sample   (if left null run in tumor_only mode)                       |
| reference_bam | Path to indexed reference BAM file                                                         |
| tumor         | Name of the tumor sample (if left null run in germline_only mode)                          |
| tumor_bam     | Path to indexed tumor BAM file                                                             |
| output_dir    | Path to the output directory. This directory will be created if it does not already exist. |
| loci          | Path to vcf file containing likely heterozygous sites (see below). Gz files supported.     |
| ref_genome_version | One of `37` or `38`. Required only when using CRAM files.                             |

The vcf file used by HMF (GermlineHetPon.37.vcf.gz) is available to download from [HMF-Pipeline-Resources](https://resources.hartwigmedicalfoundation.nl). 
The sites for hg37 were chosen by running the GATK HaplotypeCaller over 1700 germline samples and then selecting all SNP sites which are heterozygous in 800 to 900 of the samples.  The 1.3 million sites provided in this file typically result in 450k+ BAF points.  For hg38 we use GNOMAD v3 SNP sites from chr1-chrX with only a single ALT at that location and with populationAF > 0.05 and < 0.95.   This yields 7.25M sites.

Approximately 1000 sites scattered evenly through the VCF have been tagged with a SNPCHECK flag. 
The allelic frequency of these sites in the reference bam are written to the `REFERENCE.amber.snp.vcf.gz` file without any filtering to be used downstream for sample matching. 

AMBER supports both BAM and CRAM file formats. 

### Optional Arguments

| Argument              | Default | Description                                                                                       |
|-----------------------|---------|---------------------------------------------------------------------------------------------------|
| threads               | 1       | Number of threads to use                                                                          |
| min_mapping_quality   | 1       | Minimum mapping quality for an alignment to be used                                               |
| min_base_quality      | 13      | Minimum quality for a base to be considered                                                       |
| min_depth_percent     | 0.5     | Only include reference sites with read depth within min percentage of median reference read depth |
| max_depth_percent     | 1.5     | Only include reference sites with read depth within max percentage of median reference read depth |
| min_het_af_percent    | 0.4     | Minimum allelic frequency to be considered heterozygous                                           |
| max_het_af_percent    | 0.65    | Maximum allelic frequency to be considered heterozygous                                           |
| ref_genome            | NA      | Path to the reference genome fasta file. Required only when using CRAM files.                     |
| validation_stringency | STRICT  | SAM validation strategy: STRICT, SILENT, LENIENT                                                  |

### Example Usage

```
java -Xmx32G -cp amber.jar com.hartwig.hmftools.amber.AmberApplication \
   -reference COLO829R -reference_bam /run_dir/COLO829R.bam \ 
   -tumor COLO829T -tumor_bam /run_dir/COLO829T.bam \ 
   -output_dir /run_dir/amber/ \
   -threads 16 \
   -loci /path/to/GermlineHetPon.37.vcf.gz 
```

## Tumor Only Mode
If no reference BAM is supplied, AMBER will be put into tumor only mode.  In tumor only, the following behaviour is changed:
- Contamination check not run
- snp check vcf is not produced
- normalBAF and count fields are set to -1 in amber.baf.tsv
- a set of blacklisted regions (with variable germline copy number) are ignored

### Mandatory Arguments

| Argument   | Description                                                                                |
|------------|--------------------------------------------------------------------------------------------|
| tumor      | Name of the tumor sample                                                                   |
| tumor_bam  | Path to indexed tumor BAM file                                                             |
| output_dir | Path to the output directory. This directory will be created if it does not already exist. |
| loci       | Path to vcf file containing likely heterozygous sites (see below). Gz files supported.     |

### Optional Arguments

| Argument               | Default | Description                                                                   |
|------------------------|---------|-------------------------------------------------------------------------------|
| threads                | 1       | Number of threads to use                                                      |
| min_mapping_quality    | 1       | Minimum mapping quality for an alignment to be used                           |
| min_base_quality       | 13      | Minimum quality for a base to be considered                                   |
| tumor_only_min_vaf     | 0.05    | Min VAF in ref and alt in tumor only mode                                     |
| tumor_only_min_support | 2       | Min support in ref and alt in tumor only mode                                 |
| tumor_only_min_depth   | 25      | Min depth in ref and alt in tumor only mode                                   |
| ref_genome             | NA      | Path to the reference genome fasta file. Required only when using CRAM files. |

### Example Usage

```
java -Xmx32G -cp amber.jar com.hartwig.hmftools.amber.AmberApplication \
   -tumor COLO829T -tumor_bam /run_dir/COLO829T.bam \ 
   -output_dir /run_dir/amber/ \
   -threads 16 \
   -loci /path/to/GermlineHetPon.37.vcf.gz 
```
## Germline Only Mode

If the tumor / tumor bam are not specified then Amber will be run in germline only mode.   Germline mode has the following differences in behaviour
- contamination is not run
- pcf fitting is not run
- tumor fields are set to -1 for amber.baf.tsv.gz
- amber.baf.tsv.gz named with ref sample (instead of tumor sample)

## Targeted Mode

AMBER may be run on targeted data.   For more information on how to run hmftools in targeted mode please see [here](https://github.com/hartwigmedical/hmftools/blob/master/README_TARGETED.md)

## Multiple Reference / Donor mode
The `reference` and `reference_bam` arguments supports multiple arguments separated by commas. 
When run in this mode the heterozygous baf points are taken as the intersection of each of the reference bams. 
No change is made to the SNPCheck or contamination output. These will be run on the first reference bam in the list. 


## Regions of Homozygosity Algorithm
Amber outputs a file which contains continuous regions of homozygous sites.  The sex chromosomes are excluded from consideration, as are the short arms of chr 13,14,15,21 & 22 as well as regions within 1M bases of centromeric gaps and large regions of heterochromatin (ie for chr 1,chr9, chr 16).

For determination of the region each BAF site in the provided bed file is calculated as homozygous (alt or ref) or heterozygous according to the following criteria:
```
Homozygous alt:   totalReadCount==AltReadCount) OR (AltReadCount > 0.75*TotalReadCount AND POISSON.DIST(totalReadCount-AltReadCount,TotalReadCount/2,TRUE) < 0.005
Homozygous ref:  AltReadCount==0) OR (AltReadCount < 0.25*TotalReadCount AND POISSON.DIST(AltReadCount,TotalReadCount/2,TRUE) < 0.005
Heterozygous: all other sites
```
A region must meet the following criteria to be considered a region of homozygosity:
- Region must have at least 50 consecutive SNP locations with no 5 of any consecutive 50 locations heterozygous
- Region must be at least 500,000 bases in length.   For regions of <3M a soft filter of minLength is applied 
- <5% of the SNP locations in the whole region may be  heterozygous (soft filter only)

If multiple references are provided (for example donor & patient reference) then regions of homozygosity are calculated for the first reference only.

If only one chromosome is affected and the regions affected amount to more than 10Mb than mark that chromosome with UniparentalDisomy in the QC file.   Determine the ConsanguinityProportion as the sum of homozygous regions > 3MB divided by the size of the autosome genome (2.88bn bases).     Note that the expected value for consanguinityProportion for a child of direct siblings would be 0.25, dividing by a factor of 2 for each further level of removal in relationship.

## Patient Matching
The REFERENCE.amber.snp.vcf.gz contains some 1000 SNP points that can be used to identify if a new sample belongs to an existing patient. 
This is particularly important when doing cohort analysis as multiple samples from the same patient can skew results.

To enable patient matching a database is required with three tables, AmberSample, AmberMapping and AmberPatient. 
Scripts to generate these tables are available [here](../patient-db/src/main/resources/patches/amber/amber3.4_to_3.5_migration.sql).   

Each sample is loaded into AmberSample with the `LoadAmberData` application which downsamples the REFERENCE.amber.snp.vcf.gz file to 100 loci and describes each locus as:
- 1: Homozygous ref
- 2: Hetrozygous
- 3: Homozygous alt
- 0: Other (including insufficient depth (<10))

The sample is compared to all other AmberSample entries and if there is a match (>=90% of sites match), an entry is added to the AmberMapping table. 
The sample is assigned either a new patient id or an existing one (if it matches an existing sample) and the amberPatient table is updated.

A sample can be loaded with the following command:

```
java -cp amber.jar com.hartwig.hmftools.patientdb.LoadAmberData \
    -sample TUMOR \
    -amber_snp_vcf /path/to/REFERENCE.amber.snp.vcf.gz \
    -snpcheck_vcf /path/to/Amber.snpcheck.37.vcf \
    -db_user username \
    -db_pass password \
    -db_url mysql://localhost:3306/hmfpatients?serverTimezone=UTC
```

The Amber.snpcheck.37.vcf (and 38 equivalent) are available to download from [HMFTools-Resources > Amber](https://storage.googleapis.com/hmf-public/HMFtools-Resources/dna_pipeline/).

An example query to check if a sample is one of many for a patient is:
```
SELECT * FROM amberPatient WHERE patientId IN 
(SELECT patientId FROM amberPatient WHERE sampleId = 'SAMPLE');
```

The following query shows all patients with multiple samples:
```
SELECT patientId, count(*) AS sampleCount, GROUP_CONCAT(sampleId ORDER BY sampleId SEPARATOR ' ') AS samples
FROM amberPatient
GROUP BY patientId
HAVING sampleCount > 1
ORDER BY sampleCount desc;
```

## Output
| File                                 | Description                                                                              |
|--------------------------------------|------------------------------------------------------------------------------------------|
| TUMOR.amber.baf.tsv.gz               | Tab separated values (TSV) containing reference and tumor BAF at each heterozygous site. |
| TUMOR.amber.baf.pcf                  | TSV of BAF segments using PCF algorithm.                                                 |
| TUMOR.amber.qc                       | Contains QC status and comntamination rate. FAIL may indicate contamination in sample.   |
| TUMOR.amber.contamination.vcf.gz     | Entry at each homozygous site in the reference and tumor.                                |
| REFERENCE.amber.snp.vcf.gz           | Entry at each SNP location in the reference.                                             |
| REFERENCE.amber.homozygousregion.tsv | Regions of homozygosity found in the reference.                                          |

## Performance Characteristics
Performance numbers were taken from a 72 core machine using COLO829 data with an average read depth of 35 and 93 in the normal and tumor respectively. 
Elapsed time is measured in minutes. 
CPU time is minutes spent in user mode. 
Peak memory is measure in gigabytes.


| Threads | Elapsed Time | CPU Time | Peak Mem |
|---------|--------------|----------|----------|
| 1       | 144          | 230      | 15.04    |
| 8       | 22           | 164      | 18.40    |
| 16      | 12           | 164      | 21.00    |
| 32      | 8            | 170      | 21.60    |
| 48      | 7            | 199      | 21.43    |
| 64      | 6            | 221      | 21.78    |

 
# Version History and Download Links
- [3.9](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v3.9)
  - Added germline only mode support. 
  - Added `tumor_only_min_depth` argument.
  - Make `ref_genome_version` argument mandatory.
  - Added excluded SNP regions (values are hard coded) in tumor only mode.
  - Fixed a rounding issue in the classification of zygosity.
  - Removed <sample id>.amber.baf.vcf.gz output.
  - compressed baf TSV output with gzip.
- [3.8](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v3.8)
  - Added workaround for a bug in R copy_number module pcf function 
  - Fixed `NullPointerException` in tumor only mode
  - Change QC file UniparentalDisomy prints NONE instead of null
- [3.7](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v3.7)
  - Fixed `NullPointerException` when no `-ref_genome` argument is provided
- [3.6](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v3.6)
- [3.5](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v3.5)
- [3.4](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v3.4)
- [3.3](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v3.3)
- [3.2](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v3.2)
- [3.1](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v3.1)
- [3.0](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v3.0)
- [2.5](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v2.5)
- [2.4](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v2.4)
- [2.3](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v2.3)
