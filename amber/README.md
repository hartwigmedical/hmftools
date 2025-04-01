# AMBER
AMBER is designed to generate a tumor BAF file for use in PURPLE from a provided VCF of likely heterozygous SNP sites.

When using paired reference/tumor data, AMBER is also able to: 
  - detect evidence of contamination in the tumor from homozygous sites in the reference
  - facilitate sample matching / patient deduplication by recording SNPs in the germline
  - identify long regions of homozygosty and consanguinity

All resource files for this tool and the WiGiTs pipeline are available for download via the [HMF Resource page](./pipeline/README_RESOURCES.md).

The Bioconductor [copynumber](http://bioconductor.org/packages/release/bioc/html/copynumber.html) package is required for segmentation.
After installing [R](https://www.r-project.org/) or [RStudio](https://rstudio.com/), the copy number package can be added with the following R commands:
```
    library(BiocManager)
    install("copynumber")
```

AMBER requires Java 17+ to be installed.

## Paired Normal/Tumor Mode
This is the default and recommended mode.

### Mandatory Arguments

| Argument           | Description                                                                                |
|--------------------|--------------------------------------------------------------------------------------------|
| reference          | Name of the reference sample   (if left null run in tumor_only mode)                       |
| reference_bam      | Path to indexed reference BAM file                                                         |
| tumor              | Name of the tumor sample (if left null run in germline_only mode)                          |
| tumor_bam          | Path to indexed tumor BAM file                                                             |
| output_dir         | Path to the output directory. This directory will be created if it does not already exist. |
| loci               | Path to vcf file containing likely heterozygous sites (see below). Gz files supported.     |
| ref_genome_version | One of `37` or `38`                                                                        |

The loci file used by HMF for both 37 and 38 reference genomes is available to download from [HMF-Pipeline-Resources](https://resources.hartwigmedicalfoundation.nl). These loci are generated using GNOMAD v3 SNP sites (lifted over for GRCH37 version) from chr1-chrX with only a single ALT at that location and with populationAF > 0.05 and < 0.95.  These sites are further filtered to remove loci with frequently unclear zygosity in a set of 60 HMF samples, yielding around 6.3M sites overall.  

Approximately 1000 sites scattered evenly through the VCF have been tagged with a SNPCHECK flag. 
The allelic frequency of these sites in the reference bam are written to the `REFERENCE.amber.snp.vcf.gz` file without any filtering to be used downstream for sample matching. 

AMBER supports both BAM and CRAM file formats. 

### Optional Arguments

| Argument               | Default | Description                                                                                       |
|------------------------|---------|---------------------------------------------------------------------------------------------------|
| min_mapping_quality    | 50      | Minimum mapping quality for an alignment to be used                                               |
| min_base_quality       | 13      | Minimum quality for a base to be considered                                                       |
| tumor_min_depth        | 8 / 25  | Min tumor depth for a site to be considered, uses 25 in tumor-only mode                           |
| min_depth_percent      | 0.5     | Only include reference sites with read depth within min percentage of median reference read depth |
| max_depth_percent      | 1.5     | Only include reference sites with read depth within max percentage of median reference read depth |
| min_het_af_percent     | 0.4     | Minimum allelic frequency in reference sample to be considered heterozygous                       |
| max_het_af_percent     | 0.65    | Maximum allelic frequency in reference sample to be considered heterozygous                       |
| validation_stringency  | STRICT  | SAM validation strategy: STRICT, SILENT, LENIENT                                                  |
| threads                | 1       | Number of threads to use                                                                          |

### Example Usage

```
java -Xmx16G -jar amber.jar \
    -reference SAMPLE_ID_R \
    -reference_bam /sample_data/SAMPLE_ID_R.bam \ 
    -tumor SAMPLE_ID \
    -tumor_bam /sample_data/SAMPLE_ID.bam \ 
    -output_dir /sample_data/ \
    -threads 10 \
    -loci /path/to/GermlineHetPon.37.vcf.gz 
```

## Tumor Only Mode
If no reference BAM is supplied, AMBER will be put into tumor only mode.  In tumor only, the following behaviour is changed:
- Contamination check not run
- snp check vcf is not produced
- normalBAF and count fields are set to -1 in amber.baf.tsv
- a set of blacklisted regions (with variable germline copy number) are ignored
- default tumor_min_depth is raised to 25

### Mandatory Arguments

| Argument   | Description                                                                                |
|------------|--------------------------------------------------------------------------------------------|
| tumor      | Name of the tumor sample                                                                   |
| tumor_bam  | Path to indexed tumor BAM file                                                             |
| output_dir | Path to the output directory. This directory will be created if it does not already exist. |
| loci       | Path to vcf file containing likely heterozygous sites (see below). Gz files supported.     |

### Tumor-only specific optional Arguments

| Argument               | Default | Description                                   |
|------------------------|---------|-----------------------------------------------|
| tumor_only_min_vaf     | 0.05    | Min VAF in ref and alt in tumor-only mode     |
| tumor_only_min_support | 2       | Min support in ref and alt in tumor-only mode |

### Example Usage

```
java -Xmx16G -jar amber.jar \
   -tumor COLO829T -tumor_bam /run_dir/COLO829T.bam \ 
   -output_dir /run_dir/amber/ \
   -threads 16 \
   -loci /path/to/GermlineHetPon.37.vcf.gz \
   -ref_genome_version 38
```
## Germline Only Mode

If the tumor / tumor bam are not specified then Amber will be run in germline only mode.   Germline mode has the following differences in behaviour
- contamination is not run
- pcf fitting is not run
- tumor fields are set to -1 for amber.baf.tsv.gz
- amber.baf.tsv.gz named with ref sample (instead of tumor sample)

## Targeted Mode

AMBER may be run on targeted data. The differences in behaviour in Amber in targeted mode are documented here: [here](https://github.com/hartwigmedical/hmftools/blob/master/pipeline/README_TARGETED.md#amber).  

## Multiple Reference / Donor mode
The `reference` and `reference_bam` arguments supports multiple arguments separated by commas. 
When run in this mode the heterozygous baf points are taken as the intersection of each of the reference bams. 
No change is made to the SNPCheck or contamination output. These will be run on the first reference bam in the list. 

## Algorithm 

### Analysis and filtering of BAF points
When using paired reference/tumor bams, AMBER confirms these sites as heterozygous in the reference sample bam then calculates the allelic frequency of corresponding sites in the tumor bam. Only observations which meet the min map quality, min base quality and min tumor depth and with depth and AF in the specificed range are considered.    In tumor only mode, all provided sites are examined in the tumor with additional filteringon tumor vaf and allelic depth to ensure that the sites are highly unlikely to be homozygous ref or alt in the germline.
 
### Segmentation
The Bioconductor copy number package is then used to generate pcf segments from the BAF file.

### Contamination
The contamination algorithm aims to detect potential contamination from other patients’ DNA during the preparation stage.  It can only be run in Normal/Tumor mode as it works by analysing evidence in the tumor for heterozygosity at corresponding high confidence homozygous reference sites in the normal.  AMBER first detects presence of contamination and then estimates a contamination rate.

**Detect contamination**

AMBER first gathers high confidence homozygous ref sites from the normal sample – each site requires at least 7 reads supporting the REF and 0 reads supporting the ALT. The tumor sample is considered contaminated if, within these corresponding sites, 

- The number of sites with three or more ALT support reads is greater than or equals to 10 AND The number of sites with three or more ALT support reads is greater than or equals to 3% of the total number of heterozygous sites
  
OR

- The number of sites with three or more ALT support reads is greater than or equals to 10 AND The number of sites with three or more ALT support reads with VAF below 5% is greater than or equals to 0.2% of the total number of heterozygous sites

otherwise contamination is set to 0.

**Determine contamination rate**

The range of the contamination rate is evaluated from 0.001 to 1 in increments of 0.001.

In each candidate site in the Tumor sample with 2 or more ALT support
- If it is heterozygous, the expected no. contaminated reads = 0.5 * contamination rate * median read depth
- If it is homozygous with ALT, the expected no. contaminated reads = contamination rate * median read depth
- If it is homozygous with REF, the expected no. contaminated reads = 0

Assuming a Poisson distribution, the expected probability of detecting a contaminated read is
```
0.5 * Poisson(observation = ALT support, mean = ALT heterozygous mean) 
+ 0.25 * Poisson(observation = ALT support, mean = homozygous mean)
```
AMBER then excludes probabilities associated with ALT support of 0 or 1 and re-normalises the expected probability for all higher ALT support from 3 and more.

AMBER also counts the actual number of positions with ALT support. And similarly calculates the proportion belonging to each ALT support, excluding 0 and 1.

Finally, the contamination rate is chosen whose value minimises the cumulative differences between the predicted and actual proportions. 

### Regions of Homozygosity
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

### Patient Matching
The REFERENCE.amber.snp.vcf.gz contains some 1000 SNP points that can be used to identify if a new sample belongs to an existing patient. 
This is particularly important when doing cohort analysis as multiple samples from the same patient can skew results.

To enable patient matching a database is required with three tables, amberSample, amberMapping and amberPatient. 
Scripts to generate these tables are available [here](../patient-db/src/main/resources/generate_database.sql).   

Each sample is loaded into AmberSample with the `LoadAmberData` application which downsamples the REFERENCE.amber.snp.vcf.gz file to 100 loci and describes each locus as:
- 1: Homozygous ref
- 2: Hetrozygous
- 3: Homozygous alt
- 0: Other (including insufficient depth (<10))

The sample is compared to all other AmberSample entries and if there is a match (>=80% of sites match), an entry is added to the AmberMapping table. 
The sample is assigned either a new patient id or an existing one (if it matches an existing sample) and the amberPatient table is updated.

A sample can be loaded with the following command:

```
java -cp amber.jar com.hartwig.hmftools.patientdb.amber.LoadAmberData \
    -sample TUMOR \
    -amber_snp_vcf /path/to/REFERENCE.amber.snp.vcf.gz \
    -snpcheck_vcf /path/to/Amber.snpcheck.37.vcf \
    -db_user username \
    -db_pass password \
    -db_url mysql://localhost:3306/hmfpatients?serverTimezone=UTC
```

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

# Known issues / future improvements
- **Population based phasing**: Could significantly increase resolution of subclonal/low tumor fraction BAF segmentation.
 
# Version History and Download Links
- [4.0](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v4.0rc)
- [3.9](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v3.9)
- [3.8](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v3.8)
- [3.7](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v3.7)
- [3.6](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v3.6)
- [3.5](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v3.5)
- [3.4](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v3.4)
- [3.3](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v3.3)
- [3.2](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v3.2)
- [3.1](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v3.1)
- [3.0](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v3.0)
