
# PURPLE

PURPLE is a **pur**ity **pl**oidy **e**stimator. It leverages the read depth and tumor BAF to estimate the purity of a sample and generate a copy number profile.

## Required Input

### COBALT

We use COBALT to determine the number of read starts per kbase window. This must be done for both the reference and tumor samples.

For more information on how to run COBALT please refer to the [readme](https://github.com/hartwigmedical/hmftools/tree/master/count-bam-lines).


### AMBER

We use AMBER to calculate the allelic frequency in the tumor from a set of high confidence heterozygous SNPs in the reference sample.

For more information on how to run AMBER please refer to the [readme](https://github.com/hartwigmedical/hmftools/tree/master/amber).

## Usage

Argument | Default | Description
---|---|---
-ref_sample | None | Name of reference sample.
-tumor_sample | None | Name of tumor sample.
-run_dir | None | Base directory of run. Default values of output_dir, amber and cobalt will be relative to this.
-output_dir |  <run_dir>/purple | Output directory.
-baf | <amber_dir>/<tumor_sample>.amber.baf | Location of baf file. By default will look in the amber directory.
-amber | <run_dir>/amber | Location of amber directory.
-cobalt | <run_dir>/cobalt | Location of cobalt directory.
-threads | 2 | Number of threads to use.
-somatic_vcf | None | Optional location of somatic variants vcf. Sample name should match <tumor_sample>. GZ files supported.
-structural_vcf | None | Optional location of structural variants vcf. Sample name should match <tumor_sample>. GZ files supported.
-circos | None | Optional path to circos binary. When supplied, circos graphs will be written to <output_dir>/plot
-db_enabled | None | This parameter has no arguments. Optionally include if you wish to persist results to a database. Database initialization script can be found [here](https://github.com/hartwigmedical/hmftools/blob/master/patient-db/src/main/resources/generate_database.sql).
-db_user | None | Database username. Mandatory if db_enabled.
-db_pass | None | Database password. Mandatory if db_enabled.
-db_url | None | Database URL. Should be of format: `mysql://localhost:3306/hmfpatients`. Mandatory if db_enabled.

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


## Algorithm

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