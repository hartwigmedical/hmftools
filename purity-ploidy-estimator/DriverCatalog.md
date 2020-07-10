# Somatic Driver Catalog

A detailed description of our gene discovery and driver catalog is available in the [supplementary information](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-019-1689-y/MediaObjects/41586_2019_1689_MOESM1_ESM.pdf) of our ["Pan-cancer whole genome analyses of metastatic solid tumors"]((https://www.nature.com/articles/s41586-019-1689-y)) paper.

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

![Gene Classification](src/main/resources/readme/GeneClassification.png)
