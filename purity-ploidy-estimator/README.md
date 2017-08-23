
# PURPLE

PURPLE is a **pur**ity **pl**oidy **e**stimator. It leverages both CNV and BAF information to estimate the purity of a sample and then supply a copy number profile.


## R Dependencies
Copy number segmentation is done with the Bioconductor
[copynumber](http://bioconductor.org/packages/release/bioc/html/copynumber.html) package.


This can be installed in R with the following commands:
```
   source("https://bioconductor.org/biocLite.R")
   biocLite("copy number")
```


## Required Input

### Read Count

We require the number of read starts per kbase window for both the normal and tumor sample.

We use COBALT (**co**unt **ba**m **l**ines) to do this with a min quality threshold of 10.

### BAF

The BAF is the allelic frequency in the tumor of a set of high confidence heterozygous SNPs in the normal sample.

These are generated in our pipeline by the GATK haplotype caller and filtered using a bed file of 700k common germline SNPs.

In the next version, PURPLE will generate this directly from the BAMs, but for now it needs to be provided.

An example file might look:

Chromosome | Position | BAF
:--- | :---: | :---:
1 | 1186502 | 0.166
1 | 1333436 | 0.2



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