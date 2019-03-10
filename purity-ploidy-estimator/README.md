
# PURPLE

PURPLE is a **pur**ity **pl**oidy **e**stimator. It combines B-allele frequency (BAF), read depth ratios, somatic variants and structural variants to estimate the purity and copy number profile of a tumor sample.

PURPLE supports both grch 37 and 38 reference assemblies. 


## Input

The PURPLE algorithm relies on AMBER and COBALT, which are delivered with the PURPLE release, to calculate the tumor BAF and read depth ratios respectively.    

It is also strongly recommended also to run PURPLE with a high quality set of somatic SNV and INDEL calls and somatic structural variant calls.  


### COBALT

COBALT determines the read depth ratios of the supplied tumor and reference genomes. 

COBALT starts with the raw read counts per 1,000 base window for both normal and tumor samples by counting the number of alignment starts in the respective bam files with a mapping quality score of at least 10 that is neither unmapped, duplicated, secondary, nor supplementary. Windows with a GC content less than 0.2 or greater than 0.6 or with an average mappability below 0.85 are excluded from further analysis.

Next we apply a GC normalization to calculate the read ratios. We divide the read count of each window by the median read count of all windows sharing the same GC content then normalise further to the ratio of the median to mean read count of all windows. 

Finally, the reference sample ratios have a further ‘diploid’ normalization applied to them to remove megabase scale GC biases. This normalization assumes that the median ratio of each 10Mb window (minimum 1Mb readable) should be diploid for autosomes and haploid for sex chromosomes in males in the germline sample.


For more information on how to run COBALT please refer to the [readme](https://github.com/hartwigmedical/hmftools/tree/master/count-bam-lines).


### AMBER

AMBER calculates the BAF of the tumor sample by finding heterozygous locations in the reference sample from a panel of 1,344,880 common germline heterozygous SNP loci. The loci were chosen by running the GATK HaplotypeCaller over 1700 germline samples and then selecting all SNP sites which are heterozygous in 800 to 900 of the samples.

To ensure that we only capture heterozygous points, we filter the panel to only loci with allelic frequencies in the reference sample between 40% and 65% and with depth between 50% and 150% of the reference sample genome wide average. Furthermore, we filter any loci with a mapping quality < 1 or base quality < 13. This typically yields 500k-540k heterozygous germline variants per patient. 

As part of a contamination check, AMBER also finds sites in the tumor that are homologous in the reference sample using the same panel as above. A sample is considered contaminated if at least 2000 of these sites contain 3 or more reads supporting an alt in the tumor. In this case we model the expected number of non-homologous sites using a poisson distribution and estimate a contamination percent. The result of this is included in the amber QC output file.


For more information on how to run AMBER please refer to the [readme](https://github.com/hartwigmedical/hmftools/tree/master/amber).


### Structural Variant Input VCFs (optional)
Providing a high quality set of structural variant calls to PURPLE allows exact base resolution of copy number changes.   An accurate estimation of VAF at each breakend also allows PURPLE to infer copy number changes even across very short segments of the genome where a depth based estimation is inaccurate or impractical. Finally, PURPLE also supports recovery of filtered structural variant calls 

For these purposes, PURPLE provides full support and integration with the structural variant caller [GRIDSS](https://github.com/PapenfussLab/gridss). GRIDSS can be run directly on tumor and reference BAMs. Alternatively a lightweight version of GRIDSS can be used to re-analyse a set of variant calls and provide additional filtering and accurate VAF estimation.


###Somatic Variant Input VCF (optional)
An high quality set of somatic SNV and INDEL calls can also improve the accuracy and utility of PURPLE. If provided, the variants are used for enhancing the purity and ploidy fit in 2 ways.   Firstly, each solution receives a penalty for the proportion of somatic variants which have implied ploidies that are inconsistent with the minor and major allele ploidy. Secondly, for highly diploid samples, the VAFs of the somatic variants are used directly to calculate a somatic variant implied purity.

For both purposes, accurate VAF estimation is essential.   For this purpose, PURPLE requires the ‘AD’ (Alllelic Depth) field in the vcf.  High quality filtering of artefacts and false positive calls is also critical to achieving an accurate fit.


## Algorithm

There are 9 key steps in the PURPLE pipeline described in detail below:
1. Gender determination
2. Segmentation
3. Sample purity and ploidy fitting
4. Copy number smoothing
5. Inferering copy number for regions without read depth
6. Allele specific ploidy inferring
7. Recovery of structural variants and filtering of single breakends
8. Identification of germline copy number alterations that are homozygously deleted in the tumor
9. QC Status for the tumor

### 1. Gender
We examine both the AMBER and COBALT data to independently determine and validate the gender of the sample. This includes detecting the presence of Klinefelter syndrome: a chromosomal disorder resulting in 2 or more X chromosome in a male and which we have found to affect 0.2% of the male samples in our cohort. 

To determine the AMBER gender of a sample we examine the number of BAF loci outside the pseudoautosomal region of the X chromosome, anything less than 1k BAF loci is considered male. A typical female has 12-13k BAF loci on the X chromosome using our provided BED file. 

To determine the COBALT gender we first use the reference ratio to determine the number of copies of the X chromosome. A median X ratio greater than 0.65 is interpreted as 2 copies (note that nearly all female samples are very close to a ratio of 1, but a handful are significantly lower with mosaic X loss). If there is only one copy of the X chromosome the sample is male. Otherwise, we check for the presence of the Y chromosome as determined by at least 1000 data points with a median ratio > 0.05. If the Y chromosome is present (in addition to the 2 copies of the X chromosome), then the sample is male with Klinefelter syndrome. In the absence of the Y chromosome the sample is female. 

Finally we compare the AMBER and COBALT genders. If they are inconsistent we use the COBALT gender and flag the sample has having failed gender validation. 


### 2. Segmentation

We segment the genome into regions of uniform copy number by combining segments generated from the COBALT read ratios for both tumor and reference sample, the BAF points from AMBER, and passing structural variant breakpoints derived from GRIDSS. Read ratios and BAF points are segmented independently using the Bioconductor copynumber package which uses a piecewise constant fit (PCF) algorithm (with custom settings: gamma = 100, k =1). These segment breaks are then combined with the structural variants breaks according to the following rules:
1. Every structural variant break starts a new segment, as does chromosome starts, ends and centromeres. 
2. Ratio and BAF segment breaks are only included if they are at least one complete mappable read depth window away from an existing segment. 

If the segments identified by the PCF algorithm are not contiguous, then there remains some uncertainty about the actual start position of the segment. To address this, we use the PCF break as the start position but also include a min and max start position to capture the uncertainty. Segments with SV support are never uncertain. 

Once the segments have been established we map our observations to them. In each segment we take the median BAF of the tumor sample and the median read ratio of both the tumor and reference samples. We also record the number of BAF points within the segment as the BAFCount and the number of tumor read depth windows within the segment as the depth window count.

A reference sample copy number status is determined at this this stage based on the observed copy number ratio in the reference sample, either ‘DIPLOID’ (0.8<= read depth ratio<=1.2), ‘HETEROZYGOUS_DELETION’ (0.1<=ratio<0.8), ‘HOMOZYGOUS_DELETION’ (ratio<0.1),’AMPLIFICATION’(1.2<ratio<=2.2) or ‘NOISE’ (ratio>2.2). The purity fitting and smoothing steps below use only the DIPLOID germline segments.

### 3. Sample Purity and Ploidy

To estimate purity and sample ploidy, we use a model which considers a matrix of all possible sample purities and ploidies and scores each possible combination on a segment by segment basis, based on a set of principles which aim to choose the most parsimonious solution for the fit.      

The specific scoring principles applied are the following:
1. **Penalise sub-clonality**:   The major and minor allele of each segment should be close to an integer ploidy for clonal solutions.    Due to sampling noise, small deviations from integer ploidies will be observed even, but larger deviations require subclonal features and are penalised.
2. **Penalise higher ploidy solutions**:  Higher ploidies have more degenerate fits but are less biologically plausible and are given an event penalty.
3. **Penalise solutions with implausible somatic SNV ploidies**: SNVs in principle occur on only one chromatid and should not be found on both alleles.   Therefore we penalise solutions where SNV ploidies exceed the major allele ploidy.   
4. **Weigh segments by count of BAF observations**: Segments are weighted by the count of BAF observations which is treated as a proxy for confidence of BAF and read depth ratio inputs.
5. **Place more weight on segments with higher observed BAF**: segments with lower observed BAFs have more degenerate fits and are weighted less in the fit

For each [sample ploidy,purity] combination we calculate a fit score using the following formula

Fit Score = DeviationPenalty * EventPenaltyMultiplier + SomaticDeviationPenalty

The  [sample ploidy,purity] combination with the lowest fit score is selected by PURPLE as the final fit score.

Each of the 3 penalty terms is described in detail in the following sections.


#### Deviation Penalty
The deviation penalty aims to penalise [ploidy|purity] combinations which require extensive sub-clonality to explain the observed copy number pattern.

For each [ploidy|purity] combination tested an implied major and minor allele ploidy is calculated based on the observed BAF and depth ratio.    A deviation penalty is then calculated for each segment based on the implied ploidies.   The function used is designed to explicitly capture a set of intuitive rules relating to known biology of cancer genomes, specifically:
- For major allele plody > 1 and minor allele ploidy > 0 a fixed deviation penalty applies
  - the penalty depends only on the distance to the nearest integer ploidy and varies between a minimum of a small baseline deviation [0.2] and a max of 1.
  - For small deviations from an integer don’t occur any additional penalty, but once a certain noise level is exceeded the penalty grows rapidly to reflect the fact the increasing probability that the observed deviation requires a implied non-integer (subclonal) ploidy.
  - The deviation penalty is increased more slowly at lower purities reflecting the increased expected noise.   This is implemented by modeling the penalty as a normal distribution with the standard deviation a function of the purity. 
- An additional and increasing penalty multiplier applies for implied major allele ploidy < 1.   This is intended to capture the relative rarity of large homozygous deletions in cancer genomes, so potential solutions with significant amounts of homozygous deletion are penalised significantly. 
- An additional and increasing penalty also applies for implied minor allele ploidy < 0 to strongly penalise solutions which imply negative copy number and are biologically implausible.

The following chart illustrates the deviation penalty applied for each of minor and major allele ploidy at both 30% and 70% purity.

![Image Name](purity-ploidy-estimator/resources/readme/FittedPurityDeviationPenalty.png)

### Segmentation
The idea is to identify maximal set of break points so we can have discrete segments with a single BAF and absolute copy number.

1. GC normalization
2. Calculate read ratio
3. Apply piecewise constant fit (PCF) to normal ratios
4. Apply PCF to tumor ratios
5. Combine ratio segments
6. Integrate structural variant breaks (if available)

### Segment Observations
Determine median BAF and average normal and tumor depth ratios of each segment.

### Fit Purity and Ploidy

Jointly fit tumor purity and CNV ratio normalisation factor (Norm Factor) according to the following principles:

* The absolute copy number of EACH segment should be close to an integer ploidy
* The BAF of EACH segment should be close to a % implied by integer ploidy and major allele counts.
* Higher ploidies have more degenerate fits but are less biologically plausible and should be penalised
* Segments are weighted by the count of BAF observations which is treated as a proxy for segment length and confidence of BAF and ratio inputs.
* Segments with lower observed BAFs have more degenerate fits and are weighted less in the fit

An excel model of the model used to score each segment is available from [Excel Models](https://resources.hartwigmedicalfoundation.nl). 

### Absolute Copy Number Per Segment

For each segment calculate the absolute copy number implied by our fitted purity.

To allow for poly clonality BAF is ignored at this stage.

### Determine Broad CNV Regions
We want to reduce our maximal # of breakpoints into broad copy number regions across each chromosome.

At this stage we consider only regions with high BAF count and no evidence of germline noise.

Merge neighbouring regions if BAF and copy number are sufficient close to each other.

### Determine Focal Breakpoints
Extend broad regions from HC segments to the remaining LC segments where possible or create new breakpoints where necessary so as to :

* Identify potential focal amplifications/deletions between HC segments within identified broad regions
* Refine the location of breakpoints between identified broad regions
* Identify amplifications / deletions prior to the 1st HC segment or after the last HC segments


### Estimate Uncertainty
Consider all fitted purities within 10% of the best score.

Calculate min/max purity and ploidy.



## Usage

Argument | Default | Description
---|---|---
-ref_sample | None | Name of reference sample.
-tumor_sample | None | Name of tumor sample.
-run_dir | None | Base directory of run. Default values of output_dir, amber and cobalt will be relative to this.
-output_dir |  <run_dir>/purple | Output directory.
-amber | <run_dir>/amber | Location of amber directory.
-cobalt | <run_dir>/cobalt | Location of cobalt directory.
-threads | 2 | Number of threads to use.
-somatic_vcf | None | Optional location of somatic variants vcf. Sample name should match <tumor_sample>. GZ files supported.
-somatic_min_peak | 100 | Minimum number of somatic variants to consider a peak.
-somatic_min_total | 1000 | Minimum number of somatic variants required to assist highly diploid fits.
-structural_vcf | None | Optional location of structural variants vcf. Sample name should match <tumor_sample>. GZ files supported.
-sv_recovery_vcf | None | Optional location of structural variants recovery vcf. Sample name should match <tumor_sample>. GZ files supported.
-circos | None | Optional path to circos binary. When supplied, circos graphs will be written to <output_dir>/plot
-db_enabled | None | This parameter has no arguments. Optionally include if you wish to persist results to a database. Database initialization script can be found [here](https://github.com/hartwigmedical/hmftools/blob/master/patient-db/src/main/resources/generate_database.sql).
-db_user | None | Database username. Mandatory if db_enabled.
-db_pass | None | Database password. Mandatory if db_enabled.
-db_url | None | Database URL. Should be of format: `mysql://localhost:3306/hmfpatients`. Mandatory if db_enabled.
-gc_profile | None | Location of GC profile. Available to download from [HMF-Pipeline-Resources.](https://resources.hartwigmedicalfoundation.nl)
-ref_genome | Detect | Will attempt to detect reference genome from cobalt output but failing that must be either hg18 or hg38


### Example Usage

```
java -jar purple.jar \
   -run_dir "/Users/jon/hmf/analyses/COLO829" \
   -ref_sample COLO829BL \
   -tumor_sample COLO829 \
   -threads 6 \
   -gc_profile /Users/jon/hmf/analyses/COLO829/GC_profile.1000bp.cnp \
   -somatic_vcf /Users/jon/hmf/analyses/COLO829/COLO829_consensus_filtered.vcf \
   -structural_vcf /Users/jon/hmf/analyses/COLO829/somaticSV.vcf.gz \
   -circos /Users/jon/hmf/tools/circos-0.69-5/bin/circos \
   -db_enabled -db_user build -db_pass build -db_url mysql://localhost:3306/hmfpatients?serverTimezone=UTC
```
