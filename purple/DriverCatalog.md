# Somatic Driver Catalog

If a gene panel is configured, PURPLE will compile a catalog of drivers including point mutations and copy number events.  To include point mutations in the driver catalog, the somatic VCF must be annotated with SnpEff.

Running [LINX](https://github.com/hartwigmedical/hmftools/tree/master/linx) enriches the PURPLE driver catalog with the addition of both homozygous disruption drivers and fusions

## Gene Panel Configuration

The gene panel is a configuration file of which genes to add to the driver catalog and which mutation types to report in each gene.  The format is:

Field | Values | Description
---|---|---
Gene | | Name of gene
Report Missense | T/F | Report if any missense variant is found in the gene
Report Nonsense | T/F | Report if any nonsense or frameshift variant is found in the gene
Report Splice | T/F |  Report if any canonical splice acceptor or donor variant is found in the gene [+1,+2,+5,-1,-2]  Mutations affecting the last exonic base at a donor location as well as N>G variants only at the -3 acceptor base are also treated as SPLICE.
Report Amplification | T/F | Report amplification if min gene copy number > 3x sample ploidy and partial amplification if max gene copy number > 3x sample ploidy 
Report Deletion | T/F | Report if gene copy number < 0.5
Report Disruption | T/F | LINX will report ‘homozygous disruptions’ - ie disruptions where the exonic copy number of the gene is > 0.5 but where no intact copies of the gene are predicted to remain
Report Hotspot | T/F | Report somatic hotspot mutation regardless of other rules
Likelihood Type | ONCO/TSG | Calculate driver likelihood as a tumor suppressor gene or onco gene
reportGermlineVariant	| 'WILDTYPE_LOST','NONE', 'ANY','VARIANT_NOT_LOST'| Report any germline variants that meet pathogenic criteria based on specified tumor status
reportGermlineHotspot | 'WILDTYPE_LOST','NONE', 'ANY','VARIANT_NOT_LOST'| Report hotspot germline pathogenic variants based on specified tumor status


The Hartwig Medical Foundation curated gene panel is available from [HMFTools-Resources > Gene Panel](https://resources.hartwigmedicalfoundation.nl) and is updated periodically. 
A detailed description of our gene discovery and initial construction of our gene panel is available in the [supplementary information](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-019-1689-y/MediaObjects/41586_2019_1689_MOESM1_ESM.pdf) of our ["Pan-cancer whole genome analyses of metastatic solid tumors"](https://www.nature.com/articles/s41586-019-1689-y) paper.

## Gene Driver Likelihood

A driver likelihood estimate between 0 and 1 is calculated for each variant in the gene panel.
High level amplifications, deletions, fusions (when using LINX), and TERT promoter mutations are all rare so have a likelihood of 1 when found affecting a driver gene, but for coding mutations we need to account for the large number of passenger point mutations that are present throughout the genome and thus also in driver genes.

For coding mutations the driver likelihood algorithm depends on the configured 'likelihood type' in the gene panel.
The impacts of configuring `ONCO` vs `TSG` type is the following:

Characteristic | `ONCO` | `TSG`
---|---|---
High likelihood variants | Hotspot & Inframe (excluding repeat count > 7) assigned likelihood = 1 | Hotspot & Biallelic splice, indel and nonsense assigned likelihood = 1
Biallelic mutations | Ignored in dnds calculations | For biallelic mutations, biallelic TMB only used in passenger likelihood. For non biallelic, the full TMB is used
Multi-hit | Maximum likelihood used. Highest ranked variant used only in dnds calculation (missense & inframe ranked first) | Multiple mutations are additive (product of probabilities used).  For non biallelic variants highest 2 ranked variants used in dnds calculations (nonsense & splice ranked first)
Non-biallelic frameshift/Splice/Nonsense | No special treatment | driverLikelihood=max(selfLikelihood,missenseDriverLikelihood)

For point mutations, we first determine a list of variants that are highly likely to be drivers and/or highly unlikely to have occurred as passengers as driver likelihood of 1 as per the table above.
For the remaining point mutation variants these are only assigned a > 0 driver likelihood where there is a remaining excess of unallocated drivers based on the calculated dNdS rates in that gene across the cohort after applying the above rules. 
Any remaining point mutations are assigned a driver likelihood between 0 and 1 using a bayesian statistic to calculate a sample specific likelihood of each gene based on the type of variant observed (missense, nonsense, splice or INDEL) and taking into account the mutational load of the sample.   

The principle behind the likelihood method is that the likelihood of a passenger variant occurring in a particular sample should be approximately proportional to the tumor mutational burden and hence variants in samples with lower mutational burden are more likely to be drivers.

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

Note that the TMB is calculated separately for INDELs and SNVs.   For TSG dnds rates, we find that biallelic dnds rates are typically much higher than non-biallelic rates reflecting the biology that both alleles are normally required to be knocked out to cause a TSG driver.   To accommodate for this, we use count biallelic TMB in the passenger likelihood for biallelic TSG VUS.

