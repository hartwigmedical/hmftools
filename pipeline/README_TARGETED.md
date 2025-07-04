# Targeted NGS Analysis with WiGiTS

The WiGiTS tool suite fully supports targeted panel sequencing as well as WGS and WTS.  

The output from the pipeline largely matches that of WGS/WTS but with some modules excluded:
* germline analysis is typically disabled (unless the panel is also run for a normal sample)
* Cuppa is disabled
* Chord is disabled and instead Purple is used to estimate HRD
* TMB and MSI have custom calculations routines
* TPM is normalised to be in line with observed WGS rates 

For each panel a set of custom resources is required to fit and normalise the biases inherent to that specific panel. These are described in detail below.

[Oncoanalyser](https://nf-co.re/oncoanalyser/docs/usage/#panel-reference-data) can be used to run any targeted panel.

## Panel-specific resources files

To run in targeted mode, panel specific resources first need to be generated. As indicated in the below diagram, some of these are 
configured and some of these can be trained by first running `oncoanalyser` on a representative set of aligned BAMs from samples captured 
with the panel and then running a set of normalisation scripts that are described below. Here is a schematic of this process
<p align="center">
<img width="428"  alt="image" src="https://github.com/user-attachments/assets/03ec7de5-80df-4e8b-9b13-51735e37f00a">
</p>

### Manually configured files

The below files represent a basic definition of the panel and are to be created manually. Some of these files are used as inputs to the 
[panel training procedure](#training-procedure-to-build-panel-specific-resource-files). Note that the RNA resource files are only required
if your panel supports RNA data.

| Data type | File name                              | Oncoanalyser config        | Input for training? | Tool(s)  | Description                                                                                                                                                                                                                                                             |
|:----------|:---------------------------------------|:---------------------------|:--------------------|:---------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| DNA       | Driver gene panel TSV                  | `driver_gene_panel`        | Yes                 | Multiple | Defines the set of genes in the panel, and which are reported for various events (SNVs, AMPs, DELs etc). For column descriptions, see [Purple: Driver Catalog](https://github.com/hartwigmedical/hmftools/blob/master/purple/DriverCatalog.md#gene-panel-configuration) |
| DNA       | Panel regions BED                      | `target_region_bed`        | Yes                 | Multiple | Standard BED file defining the panel regions. Other regions are ignored in the training. May include gene exons or other regions of interest. The bed file must be sorted and exclude any ALT contigs or non standard human chromosomes                                 |
| DNA       | MSI indels TSV                         | `target_region_msi_indels` | Yes                 | PURPLE   | List of MSI loci to consider in MSI model. Columns: Chromosome, Position                                                                                                                                                                                                |
| DNA       | TMB/MSI configuration TSV <sup>1</sup> | `target_region_ratios`     |                     | PURPLE   | Configuration for normalising TMB/MSI values from panel levels to WGS levels                                                                                                                                                                                            |
| RNA       | RNA panel genes CSV                    | `isofox_gene_ids`          | Yes                 | ISOFOX   | Ensembl gene IDs and gene names for genes in the RNA panel. Note: these genes may not necessarily match genes covered by `target_region_bed`                                                                                                                            |
| RNA       | Expected counts CSV                    | `isofox_counts`            |                     | ISOFOX   | Pre-computed expected counts per transcript and gene. Recommended to use the `read_151_exp_counts.*.csv` from the [WiGiTS reference data](https://nf-co.re/oncoanalyser/docs/usage/#reference-data-urls)                                                                |
| RNA       | Expected GC ratios CSV                 | `isofox_gc_ratios`         |                     | ISOFOX   | Pre-computed expected GC ratio counts per transcript. Recommended to use the `read_100_exp_gc_ratios.*.csv` from the [WiGiTS reference data](https://nf-co.re/oncoanalyser/docs/usage/#reference-data-urls)                                                             |

1. The default TMB/MSI configuration of this file is shown below:
```
TmbRatio	TmlRatio	MsiIndelRatio	Msi23BaseAF	Msi4BaseAF	CodingBaseFactor
0.05	0.74	220	0.15	0.08	150000
```

### Output files training procedure

The following files are the output of the training process described below:

| Data type | File name                        | Config name                   | Tool(s) | Description                                                                                                        |
|:----------|:---------------------------------|:------------------------------|:--------|:-------------------------------------------------------------------------------------------------------------------|
| DNA       | Target regions normalisation TSV | `target_region_normalisation` | COBALT  | Normalise copy number regions and mask off-target regions                                                          |
| DNA       | Panel artefact PON TSV           | `pon_artefacts`               | PAVE    | Panel-specific PON to apply in addition to standard WGS PON                                                        |
| RNA       | TPM normalisation CSV            | `isofox_tpm_norm`             | ISOFOX  | Normalise TPM for genes in panel (only required for panels with targeted RNA content). Not required for PANEL+WTS) |

These files are then used by the pipeline in panel mode to produce well-calibrated, accurate results ideally without panel-specific biases.

## Training procedure to build panel-specific resource files

An initial set of input samples, recommended to number at least 20, are use to 'train' the pipeline. 
In particular this identifies variance in read depth and variant calling compared with whole genome and transcriptome. 
Note that the assumption of this training is that the _median_ relative copy number for any given gene should not deviate too systematically from the ploidy across the cohort. 
In general this is a relatively safe assumption for a pan-cancer dataset, but in cancer specific training sets, there may be certain recurrently copy-number-altered genes that violate this assumption.   

There are 2 steps in the training procedure:

### STEP 1: Run COBALT, AMBER, SAGE & ISOFOX on the targeted samples in WGTS mode

This can be done by running the samples through oncoanalyser in WGTS mode.

Alternatively COBALT, AMBER & SAGE & ISOFOX can be run on each sample with the below configurations

COBALT
```
java -jar cobalt.jar \
    -tumor SAMPLE_ID \
    -tumor_bam /sample_data/SAMPLE_ID.bam \ 
    -output_dir /sample_data/cobalt \ 
    -gc_profile /ref_data/GC_profile.1000bp.37.cnp \
    -tumor_only_diploid_bed DiploidRegions.37.bed.gz \
    -threads 10 \ 
```
AMBER
Amber is required to determined the gender of the samples for copy number normalisation.  Note that we recommend to lower the tumor_min_depth to ensure sufficient coverage of heterozygous points.
```
java -jar amber.jar com.hartwig.hmftools.amber.AmberApplication \
    -tumor SAMPLE_ID \
    -tumor_bam /sample_data/SAMPLE_ID.bam \
    -loci /path/to/GermlineHetPon.37.vcf.gz \ 
    -tumor_min_depth 2 \ 
    -output_dir /sample_data/ \
    -threads 10 \
```
SAGE
```
java -jar sage.jar \
    -tumor SAMPLE_ID \ 
    -tumor_bam /sample_data/SAMPLE_ID.bam \
    -ref_genome_version 37 \
    -ref_genome /ref_data/refGenome.fasta \
    -hotspots /ref_data/KnownHotspots.37.vcf.gz \
    -panel_bed /path/to/ActionableCodingPanel.37.bed.gz \
    -high_confidence_bed /ref_data/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed \
    -ensembl_data_dir /ref_data/ensembl_cache/ \
    -output_vcf /sample_data/COLO829v003.sage.vcf.gz \
    -threads 16 \
```

ISOFOX (panels with targeted RNA only)
```
java -jar isofox.jar \
    -sample SAMPLE_ID \
    -functions TRANSCRIPT_COUNTS \
    -bam_file /sample_data/SAMPLE_ID.RNA.bam \ 
    -ref_genome /ref_data/ref-genome.fasta \
    -ensembl_data_dir /ref_data/ensembl_cache/ \ 
    -exp_counts_file /ref_data/read_151_exp_counts.csv \ 
    -exp_gc_ratios_file /ref_data/read_100_exp_gc_ratios.csv \ 
    -output_dir /sample_data/ \
    -threads 16 \
```

Note: Isofox requires the expected counts file to have been generated for the correct RNA read length. See Isofox read-me for details.

### STEP 2: Run training normalisation scripts

#### Target Regions Normalisation TSV

Run the Cobalt normalisation file builder command described below on the training samples output.  This performs the following steps
- for each 1K region covering any target region, extract each sample's tumor read count and the GC profile mappability and GC ratio bucket
- calculate median and median read counts for each sample, and overall sample mean and median counts
- normalise each sample's tumor read counts per region
- calculate a median read count from all samples per GC ratio bucket
- write a relative enrichment for each region to the output file, with a min enrichment of 0.1
- if no WGS is available for normalisation, the tumorGCRatio is assumed to be 1 for autosomes. The gender of each sample must be provided. Female samples are excluded from Y chromosome normalisation and males use a tumorGCRatio of 0.5 for the sex chromosomes

The output of this process is a target regions normalisation file with the expected relative enrichment for each on target region.

##### Arguments
Field | Description
---|---
sample_id_file | CSV with SampleId column header
cobalt_dir | Cobalt output directory from step 1 above
amber_dir | Amber output directory from step 2 above
ref_genome_version | V37 or V38
gc_profile | As used in Cobalt and Purple
target_regions_bed | Definition of target regions
output_file | Output normalisation TSV file

##### Command
```
java -cp cobalt.jar com.hartwig.hmftools.cobalt.norm.NormalisationFileBuilder 
  -sample_id_file sample_ids.csv \
  -cobalt_dir /training_sample_data/ \
  -amber_dir /training_sample_data/ \
  -ref_genome_version V37 \
  -gc_profile /ref_data/GC_profile.1000bp.37.cnp \
  -target_regions_bed /ref_data/target_regions_definition.37.bed \
  -output_file /ref_data/target_regions.cobalt_normalisation.37.tsv \
  -log_debug \
```

#### Panel Artefact PON
In addition to the WGS PON, Pave utilises a panel-specific PON to capture panel specific artefacts. 
Any non-hotspot variant found 3 or more times with a qual (TQP) of > 40 and modified map-qual factor of > -10 is added to the panel PON.

To generate this file all the Pave PonBuilder to make the additional PON file:

```
java -cp pave.jar com.hartwig.hmftools.pave.pon_gen.PonBuilder \
  -sample_id_file training_sample_ids.csv \
  -vcf_path "/training_sample_data/*.sage.vcf.gz" \
  -ref_genome_version V38 \
  -output_dir /ref_data/ \
  -log_debug \
```

#### TPM normalisation

The TPM normalisation is only required for panels with RNA coverage (eg. TSO500) in order to normalise the TPM so that it is equivalent to WTS.  This can be generated using the Isofox Normalisation Builder

```
java -cp isofox.jar com.hartwig.hmftools.isofox.cohort.CohortAnalyser \
  -root_data_dir /sample_data/ \
  -sample_data_file /sample_data/training_sample_ids.csv \
  -analyses PANEL_TPM_NORMALISATION \
  -gene_distribution_file /ref_data/isofox.wgs_gene_distribution.37.csv \
  -gene_id_file /ref_data/targeted_pane_gene_ids.csv \
  -output_dir /ref_data/ \
```

The median adjusted TPM for each gene across the panel samples. If the median is zero, a replacement value of 0.01 is used instead. The adjustment factor is calculated by dividing the panel median value by the corresponding whole genome value for each gene.

Note: The adjustment factors are calculated at the gene level and not at the transcript level. This means the adjusted TPMs for transcripts from panel sequencing are not reliable.

## Pipeline Tool Functional Differences

<TO DO: make into a table>

### Cobalt
Cobalt normalises copy number and masks off-target regions according to the CN normalisation file. If a targetRegions file is provided, then a target enrichment rate is calculated simply as the median tumorGCRatio for the specified regions.
Any depth windows outside of the targetRegions file are masked so that they are ignored downstream by PURPLE. Depth windows found in the TSV file are normalised first by the overall target enrichment rate for the sample, then by the relativeEnrichment for that depth window and finally by the normal GC bias adjustment. The GC bias is calculated using on target regions only.

### Amber
The following filters are applied:
* min_depth (in tumor) > 25
* Tumor ref and alt support >= 2
* Min_depth_percent and max_depth_percent are not applied
* Tumor ref and alt VAF >= 0.05
* Amber loci must be within 300 bases of a target region

### Purple
To estimate MSI, a set of microsatellites with high coverage in the panel must also be defined.

#### MSI estimate
For a set of microsatellite sites defined in the MSI target bed file count the number of passing variants at MSI sites ignoring SNV, MNV and 1 base deletes and requiring a VAF cutoff of > 0.15 for 2 and 3 base deletes or 0.08 for 4+ base deletes or any length insertion.

We estimate MSI rate as:
```
MSIndelsPerMb = 220 * # of MSI variants / # of MSI sites in panel
```

#### TML & TMB estimate
A custom model is used for TMB estimated in targeted mode. The main challenges of the model is to determine variants are included in the TMB estimate. PURPLE selects variants that meet the following criteria:
- Coding effect <> NONE
- GNDFreq <0.00005
- GENE in PANEL and not in {HLA-A,HLA-B,HLA-C,PIM1,BCL2}
- Type = SNV
- !HOTSPOT
- 0.05 <AF < 0.9

Each variant included is classified as ‘somatic’ if somatic likelihood = HIGH.    If somatic likelihood = MEDIUM, then the variant is marked as 'unclear'.

The final somatic count estimate is set to = somatic + unclear^2 / ( CodingBases/170,000 + unclear).

This function is intended to reflect that when the number of unclear variants is less than expected germline variants then most unclear variants will be germline, whereas where the number of unclear variants is very high most will be somatic.

Using this number we then estimate the mutational burden as follows
```
TML = somatic Variant Estimate / CodingBases * RefGenomeCodingBases 
TMB = 0.05 * TML + MSIIndelPerMb
```
The 0.05 conversion from TML to TMB is the empirically observed relationship in the Hartwig database.

For driver likelihood calculations, we assume 20% of variants are biallelic for targeted sequencing samples.

#### Purple

The following special rules apply to the consrtuction of the driver catalog
- **DELS**: Don’t report DELS >10Mb or if the copy number segment has less than 3 depth windows (unless supported by SV on both sides)
- **PARTIAL_AMP**: only in genes with known pathogenic exon deletions {BRAF, EGFR, CTNNB1, CBL,MET, ALK, PDGFRA}

There is also no somatic fit mode or somatic penalty and no SV recovery in PURPLE in targeted mode.

### Isofox
TPM is normalised to bring panel gene expression in-line with WGS expression rates.

### Sage
Sage is run in high depth mode. See Sage readme for details.

## Recommended parameter values
The following parameters are calibrated for panel sequencing and are set differently to WGS. These are the default panel parameter values.

Amber

```
-target_regions_bed Target Regions BED file
```

Cobalt
```
-target_region Target Regions Normalisation TSV
-pcf_gamma 50
```

Sage
```
-high_depth_mode
```

Purple
```
-target_regions_bed Target Regions BED file
-target_regions_ratios Target Regions Ratios file
-target_regions_msi_indels Target Regions MSI Indels file
```

## Future improvements
* **HRD** - vCHORD (HRD classifier) is completed and will be made available in an upcoming release.
* **TOO** - vCUPPA (CNN based tissue of origin classifier) is under development
