# Somatic Driver Catalog

A detailed description of our gene discovery and driver catalog is available in the [supplementary information](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-019-1689-y/MediaObjects/41586_2019_1689_MOESM1_ESM.pdf) of our ["Pan-cancer whole genome analyses of metastatic solid tumors"](https://www.nature.com/articles/s41586-019-1689-y) paper.

Since publication we have have made a number of significant improvements including:
- Dnds and significant amps and dels updated to run on larger cohort of 3524 samples
- Total number of driver genes extended from 380 to 408 genes
- Reportable point mutations extended to include all genes with actionable point mutations (B evidence level or higher)
- Reported amps and dels extended to include all actionable amps and dels (B evidence level or higher)
- TSG / Oncogene classifications updated to mark genes dominated by INDELs as TSG, and to use known and actionable amplifications and deletions to assist classification for COSMIC ambiguous, non significant genes.
- Biallelic TSG missense variants no longer reported as driverlikelihood = 1.  Instead dnds is calculated separately for bialllelic and non-biallelic variants.

These are described in more detail below. 

## PURPLE Driver Catalog

PURPLE outputs a driver catalog that is automatically populated with significant amplifications (mininum exonic copy number > 3 * sample ploidy) and deletions (minimum exonic copy number < 0.5). 
Amplifications are restricted to the target amplification genes described in the supplementary information and any oncogenic point mutation genes.
Deletions are restricted to the target deletion genes described in the supplementary information and any tumor suppressor point mutation genes.

To include point mutations in the driver catalog, the somatic VCF must be annotated with SnpEff. 

The PURPLE driver catalog does NOT include fusions. 
Running [LINX](https://github.com/hartwigmedical/hmftools/tree/master/sv-linx) can improve the PURPLE driver catalog with the addition of fusion drivers. 

# Improvements

The following sections highlight the changes from what is described in the supplementary information. 

## Significantly amplified & deleted driver gene discovery changes

We updated the method described in the supplementary information to account for the larger cohort by raising the cutoffs for significant deletions to >7 homozygous deletes (previously >5) and for amplifications to a score of >35 (previously >29).
For annotation of unknown ambiguous amplification and deletion targets, we modified the prioritisation so that the longest protein coding gene was selected as the target gene ahead of the highest scoring gene.     
Target regions with more than 50% of the amplifications or deletions bounded by a centromere or telomere are also marked as telomeric or centromeric.

## Panel of driver genes for point mutations

We created a gene panel for reporting point mutations using the union of
- Martincorena significantly mutated genes (filtered to significance of q<0.01)
- HMF significantly mutated genes (q<0.01) at a global level or at cancer type level
- Cosmic Curated Genes (v83)
- Genes with actionable variants or ranges at B evidence level or higher from our knowledgebase tool (using the UNION of CIVIC, CGI & OncoKb).  

## Determine TSG or Oncogene status of each gene in point mutation panel

With the addition of new genes in the panel we reran our gene classification logic.

We used a logistic regression model to classify the genes in our pane as either tumor suppressor gene (TSG) or oncogene. 
We trained the model using unambiguous classifications from the Comic census genes, i.e. a gene was considered either an Oncogene or TSG but not both. 
We determined that the dNdS missense and nonsense ratios (w_missense and w_nonsense) are both significant predictors of the classification. 
The coefficients are given in the table below. 

 Type | Estimate | Std. Error | z value | Pr(>z)   
---|---|---|---|---
intercept | -0.0995 | 0.3525 | -0.282 | 0.7778 
w_missense | -0.5644 | 0.2282 -2.473 | 0.0135
w_nonsense | 0.5648 | 0.1198 | 4.714 | 2.43e-06


We applied the model to all significantly mutated genes in HMF as well as any ambiguous Cosmic census genes with the following exceptions:
- Ambiguous COSMIC census genes which were not significant in HMF but which have known or actionable amplifications were marked as Oncogenes
- Ambiguous COSMIC census genes which were not significant in HMF but which have known or actionable deletions were marked as TSG.
- Ambiguous COSMIC census genes which were not significant in HMF but which had >50% of calculated drivers as INDELs were marked as TSG.

The following figure shows all genes that have classified using the logistic regression model. 
Figures A and C show the likelihood of a gene being classified as a TSG under a single variate logistic model of w_missense and w_nonsense respectively. 
Figure B shows the classification after the multivariate regression using both predictors. 

<p align="center">
    <img src="src/main/resources/readme/GeneClassification.png" width="600" alt="Gene Classification">
</p>

## Gene Driver Likelihood

A driver likelihood estimate between 0 and 1 is calculated for each variant in the gene panel. 
High level amplifications, Deletions, Fusions (when using LINX), and TERT promoter mutations are all rare so have a likelihood of 1 when found affecting a driver gene, but for coding mutations we need to account for the large number of passenger point mutations that are present throughout the genome and thus also in driver genes.

For coding mutations we also mark coding mutations that are highly likely to be drivers and/or highly unlikely to have occurred as passengers as driver likelihood of 1, specifically:
- Known hotspot variants
- Inframe indels in oncogenes with repeat count < 8 repeats. Longer repeat count contexts are excluded as these are often mutated by chance in MSI samples
- Biallelic splice, nonsense or indel variants in tumor suppressor genes

For the remaining point mutation variants (non-hotspot missense variants in oncogenes and non-hotspot variants of all types in TSG) these are only assigned a > 0 driver likelihood where there is a remaining excess of unallocated drivers based on the calculated dNdS rates in that gene across the cohort after applying the above rules. 
Any remaining point mutations are assigned a driver likelihood between 0 and 1 using a bayesian statistic to calculate a sample specific likelihood of each gene based on the type of variant observed (missense, nonsense, splice or INDEL) and taking into account the mutational load of the sample.
For TSG dnds rates, we find that biallelic dnds rates are typically much higher than non-biallelic rates reflecting the biology that both alleles are normally required to be knocked out to cause a TSG driver.   
Hence for missense variants in TSG where we observed a higher biallelic dnds rate compared to non-biallelic, we treat bialllelic and non-biallelic missense variants as independent variant classes and calculated dnds rates, passenger and driver counts separately. 
If the biallelic dnds rate is lower than the non-biallelic rate for a TSG then all missense variants for that gene are pooled together and the combined dnds rate for that gene is used.

The principle behind the likelihood method is that the likelihood of a passenger variant occuring in a particular sample should be approximately proportional to the tumor mutational burden and hence variants in samples with lower mutational burden are more likely to be drivers.   
The sample specific likelihood of a residual excess variant being a driver is estimated for each gene using the following formula:

```
P(Driver|Variant) = P(Driver) / (P(Driver) + P(Variant|Non-Driver) * (1-P(Driver)))
```

where P(Driver) in a given gene is assumed to be equal across all samples in the cohort, ie:

```
P(Driver) = (residual unallocated drivers in gene) / # of samples in cohort
```

And P(Variant|Non-Driver), the probability of observing n or more passenger variants of a particular variant type in a sample in a given gene, is assumed to vary according to tumor mutational burden, and is modelled as a poisson process:

```
P(Variant|Non-Driver) = 1 - poisson(λ = TMB(Sample) / TMB(Cohort) * (# of passenger variants in cohort),k=n-1)
```

Note that the TMB for both the samples and cohort are calculated separately for INDELs and SNVs, and in the case of TSG separately for bialllelic and non biallelic variants.

