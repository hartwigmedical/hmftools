**WORK IN PROGRESS**

# Targeted NGS Analysis in HMF Tools

TO DO:    Aims to generically handling targeted sequencing

Handles TSO500 + HMF
link to demo


## Special resources files for targeted mode

File Name | Tool | Purpose
-- | -- | --

TO DO

## Recommended configuration changes in targeted mode

SAGE
```
-hotspot_min_tumor_vaf 0.01
-hotspot_min_tumor_qual 150
-panel_min_tumor_qual 250
-high_confidence_min_tumor_qual 350
-low_confidence_min_tumor_qual 500
```
GRIPSS
```
-hard_min_tumor_qual 200 
-min_qual_break_point 1000
-min_qual_break_end 1000
-filter_sgls
```

## Targeted specific methods
### Cobalt depth coverage normalisation
#### Generation of targetRegions CN normalisation file
A tsv file used for COBALT targeted CN normalisation can be prepared from a set of matching tumor targeted BAMs, + optionally a WGS BAM (if available) and a target regions bed file. The process for generating the tsv file is as follows:

1. Set targetRegions = cobalt 1kb depth windows that meet normal GC and mappability criteria and overlap or partially overlap at least 1 region specified in the primary targets bed file

2. Run cobalt on bams from targeted and matching WGS samples. For the targeted samples, calculate the targeted regions enrichment rate as median(tumorGCRatio) of the targetRegions.

For each depth window in targetRegions, calculate relativeEnrichment across the n samples as:
```
relativeEnrichment = median[tumorGCRatio(Targeted(i))/tumorGRatio(WGS(i)) / targetEnrichmenTRate(i)] 

If relativeEnrichment < 0.1 then the region is masked.
```
If no WGS is available for normalisation, the tumorGCRatio is assumed to be 1 for autosomes. The gender of each sample must be provided. Female samples are excluded from Y chromosome normalisation and males use a tumorGCRatio of 0.5 for the sex chromosomes

3. Write targetRegions normalisation file:
- chromosome
- positonStart
- relativeEnrichment

#### Calculation of target enrichment rate when targetRegions tsv specified
If a targetRegions file is provided, then a target enrichment rate is calculated simply as the median tumorGCRatio for the specified regions
Masking and normalisation of GC ratio when targetRegions tsv specified
If a targetRegions tsv file is provided then any depth windows outside of the targetRegions file are masked so that they are ignored downstream by PURPLE. Depth windows found in the tsv file are normalised first by the overall target enrichment rate for the sample and then by the relativeEnrichment for that depth window.

#### Off target normalisation
For each 100kb bucket that does not overlap an on-target region and has at least half of the depth windows that meet the COBALT GC and mappability criteria, calculate the median tumor ratios. Normalise such that the median of all 100kb buckets for the bam is 1. Note we don’t use the mean since some regions still contain dna that has been highly enriched by the targeted panel which could skew the average.

Note, the off target normalisation is not currently used in the targeted output.

### PURPLE MSI 

For a set of microsatellite sites defined in the MSI target bed file count the number of passing variants at MSI sites ignoring SNV, MNV and 1 base deletes and requiring a VAF cutoff of > 0.15 for 2 and 3 base deletes or 0.08 for 4+ base deletes or any length insertion.

We estimate MSI rate as:
```
MSIndelsPerMb = 220 * # of MSI variants / # of MSI sites in panel
```

### PURPLE TML & TMB estimate

A custom model is used for TMB estimated in targeted mode. The main challenges of the model is to Variants are included in the TMB estimate that meet the following criteria:
- Coding effect <> NONE
- GNDFreq <0.00005
- GENE in PANEL and not in {HLA-A,HLA-B,HLA-C,PIM1,BCL2} 
- Type = SNV
- !HOTSPOT
- AF < 0.9

Each variant included is classified as ‘somatic’ if AF is more than 0.08 away from the expected germline allele frequencies based on the minor and major allele copy numbers at that loci and the purity of the sample. Other variants are classified as ‘unclear’

The final somatic count estimate is set to = somatic + unclear^2 / ( CodingBases/170,000 + unclear).

This function is intended to reflect that when the number of unclear variants is less than expected germline variants then most unclear variants will be germline, whereas where the number of unclear variants is very high most will be somatic.

Using this number we then estimate the mutational burden as follows
```
TML = 0.74 * somatic Variant Estimate / CodingBases * RefGenomeCodingBases 
TMB = 0.05 * TML + MSIIndelPerMb
```
The constant 0.74 is the approximate proportion of coding variants expected to be missense, since TML is defined as count of missense variants in Hartwig’s platform. The 0.05 conversion from TML to TMB is the empirically observed relationship in the Hartwig database.


### PURPLE driver catalog differences

TO DO

## Panel specific PONs	

We constructed a panel specific PON separately for each panel based on the 48 HMF and XX TSO500 samples we have analysed. Any non hotspot variant found 3 or more times with a qual of > 100 was excluded.

We also added the following HOTSPOT variants to both PONs:
chr3:142555897:AT>A
chrX:67545400:GGCA>G

Variants that were already in our normal PON were excluded.

## Future improvements

- **Overlapping reads** - Overlapping reads in the same fragment are common in FFPE Panel libraries. We should ideally only use the highest qual at each base in case of overlap.
- **UMI** - We should switch to using UMI to better filter duplicates
- **Off Target normalisation and integration** - This is implemented, but not used as currently does not yield a benefit over on target alone.
MSI thresholds - We could better estimate if we had a more diverse range of samples for testing with known MSIndelsPerMb rates around and above the MSI cutoff.
- **HRD prediction** - We can likely train a custom model, but we don’t have enough known HRD+ samples with panel data at the moment to evaluate.
- **Purity & ploidy estimates** - Purity and ploidy estimates are only correct approximately half the time.   The fit could be made more robust by improving -COBALT/AMBER parameterisation, merging on and off target regions or changing PURPLE fit functionality
