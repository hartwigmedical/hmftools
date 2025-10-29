# Targeted NGS analysis with WiGiTS

In addition to WGS/WTS data, the WiGiTS tool suite fully supports analysis of targeted panel sequencing data via
[Oncoanalyser](https://nf-co.re/oncoanalyser/docs/usage). This document provides instructions for creating the 
[panel-specific resource files](#panel-specific-resource-files) that are required before `oncoanalyser` can be run on your custom panel 
samples.

The targeted pipeline largely matches the WGS/WTS pipeline but with some modules excluded:
* Germline analysis is typically disabled (unless the panel is also run for a normal sample)
* [CUPPA](https://github.com/hartwigmedical/hmftools/tree/master/cuppa) is disabled (for tissue of origin prediction)
* [CHORD](https://github.com/hartwigmedical/hmftools/tree/master/chord) is disabled and instead [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/chord) is used to estimate HRD
* TMB and MSI have custom calculations routines
* TPM is normalised to be in line with observed WGS rates

## Table of contents
<!-- TOC -->
* [Panel-specific resource files](#panel-specific-resource-files)
  * [Manually configured files](#manually-configured-files)
  * [Output files from training procedure](#output-files-from-training-procedure)
* [Panel resources training procedure](#panel-resources-training-procedure)
  * [Running the training procedure](#running-the-training-procedure)
  * [Overview of training procedure](#overview-of-training-procedure)
    * [Step 1: Run COBALT, AMBER, SAGE & ISOFOX](#step-1-run-cobalt-amber-sage--isofox)
    * [Step 2: Run training normalisation scripts](#step-2-run-training-normalisation-scripts)
* [Targeted analysis pipeline differences](#targeted-analysis-pipeline-differences-)
  * [Functional differences per tool](#functional-differences-per-tool)
  * [Recommended parameter values per tool](#recommended-parameter-values-per-tool)
* [Future improvements](#future-improvements)
<!-- TOC -->

## Panel-specific resource files

The [panel-specific resource files](#panel-specific-resource-files) fits and normalises the biases inherent to your panel. 
[Manual configuration](#manually-configured-files) is required for some of these files, which in turn are used as input for a 
[training procedure](#panel-resources-training-procedure) that generates the remaining panel-specific resource files 
from a representative set of sample BAMs or FASTQs from your panel (**≥20 samples recommended**). The below diagram summarises the 
generation of panel-specific resource files.

<p align="center">
<img width="428"  alt="image" src="https://github.com/user-attachments/assets/03ec7de5-80df-4e8b-9b13-51735e37f00a">
</p>

### Manually configured files

The below files represent a basic definition of the panel and are to be created manually. All files 
except for the TMB/MSI configuration TSV are used as inputs to the [panel training procedure](#panel-resources-training-procedure).

| Data type | File name                                              | Oncoanalyser config        | Input for training? | Tool(s)  | Description                                                                                                                                                                                                                                                             |
|:----------|:-------------------------------------------------------|:---------------------------|:--------------------|:---------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| DNA       | [Driver gene panel TSV](#driver-gene-panel-tsv)        | `driver_gene_panel`        | Yes                 | Multiple | Defines the set of genes in the panel, and which are reported for various events (SNVs, AMPs, DELs etc). For column descriptions, see [Purple: Driver Catalog](https://github.com/hartwigmedical/hmftools/blob/master/purple/DriverCatalog.md#gene-panel-configuration) |
| DNA       | [Panel regions BED](#panel-regions-bed)                | `target_region_bed`        | Yes                 | Multiple | Standard BED file defining the panel regions. Other regions are ignored in the training. May include gene exons or other regions of interest. The bed file must be sorted and exclude any ALT contigs or non standard human chromosomes                                 |
| DNA       | [MSI indels TSV](#msi-indels-tsv)                      | `target_region_msi_indels` | Yes                 | PURPLE   | Chromosome/positions of MSI loci to consider in MSI model                                                                                                                                                                                                               |
| DNA       | [TMB/MSI configuration TSV](#tmbmsi-configuration-tsv) | `target_region_ratios`     |                     | PURPLE   | Configuration for normalising TMB/MSI values from panel levels to WGS levels. See below for recommended defaults                                                                                                                                                        |
| RNA       | [RNA panel genes CSV](#rna-panel-genes-csv)            | `isofox_gene_ids`          | Yes                 | ISOFOX   | Ensembl gene IDs and gene names for genes in the RNA panel. Note: these genes may not necessarily match genes covered by `target_region_bed`                                                                                                                            |

Below are example snippets of each manually configured file. For full examples of these files, please see the 
[TSO500 panel resource files](https://nf-co.re/oncoanalyser/docs/usage/#reference-data-urls) from `oncoanalyser`.  

#### Driver gene panel TSV
```
gene	reportMissense	reportNonsense	reportSplice	reportDeletion	reportDisruption	reportAmplification	reportSomaticHotspot	likelihoodType	reportGermlineVariant	reportGermlineHotspot	reportGermlineDisruption	additionalReportedTranscripts
CDKN2A	TRUE	TRUE	TRUE	TRUE	TRUE	FALSE	TRUE	TSG	NONE	NONE	FALSE	ENST00000579755
CDKN2C	TRUE	TRUE	TRUE	TRUE	TRUE	FALSE	TRUE	TSG	NONE	NONE	FALSE	
...
```

#### Panel regions BED
```
1	2488103	2488172	TNFRSF14_1_CODING
1	2489164	2489273	TNFRSF14_2_CODING
...
```

#### MSI indels TSV

This file have empty rows (i.e. only the header) if you do not expect any of your samples to have MSI.

```
Chromosome	Position
1	16200729
...
```

#### TMB/MSI configuration TSV

This example shows the recommended defaults that can be used if you are unsure of how to set these values.

```
TmbRatio	TmlRatio	MsiIndelRatio	Msi23BaseAF	Msi4BaseAF	CodingBaseFactor
0.05	0.74	220	0.15	0.08	150000
```

#### RNA panel genes CSV
```
GeneId,GeneName
ENSG00000097007,ABL1
...
```

### Output files from training procedure

The following files are the output of the training process described below:

| Data type | File name                        | Config name                   | Tool(s) | Description                                                                                                        |
|:----------|:---------------------------------|:------------------------------|:--------|:-------------------------------------------------------------------------------------------------------------------|
| DNA       | Target regions normalisation TSV | `target_region_normalisation` | COBALT  | Normalise copy number regions and mask off-target regions                                                          |
| DNA       | Panel artefact PON TSV           | `pon_artefacts`               | PAVE    | Panel-specific PON to apply in addition to standard WGS PON                                                        |
| RNA       | TPM normalisation CSV            | `isofox_tpm_norm`             | ISOFOX  | Normalise TPM for genes in panel (only required for panels with targeted RNA content). Not required for PANEL+WTS) |

These files are then used by the pipeline in panel mode to produce well-calibrated, accurate results ideally without panel-specific biases.

## Panel resources training procedure

The panel training procedure identifies variance in read depth and variant calling compared with whole genome and transcriptome.

> [!NOTE]
> An assumption of the training procedure is that the _median_ relative copy number for any given gene should not deviate too systematically 
> from the ploidy across the cohort. In general, this is a relatively safe assumption for a pan-cancer dataset, but in cancer specific 
> training sets, there may be certain genes with recurrent copy number alterations that violate this assumption.

### Running the training procedure
The training procedure can be run with `oncoanalyser` using `--mode panel_resource_creation`, First create a samplesheet with a set of 
representative samples (**≥20 samples recommended**) from your panel sequencing run:

```
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.bam
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bai,/path/to/PATIENT1-T.dna.bam.bai
PATIENT2,PATIENT2,PATIENT2-T,tumor,dna,bam,/path/to/PATIENT2-T.dna.bam
PATIENT2,PATIENT2,PATIENT2-T,tumor,dna,bai,/path/to/PATIENT2-T.dna.bam.bai
```

Then, run `oncoanalyser` with `--mode panel_resource_creation` (`--isofox_*` arguments only required if panel supports RNA data):

```bash
nextflow run nf-core/oncoanalyser \
  -profile docker \
  -revision 2.2.0 \
  -config refdata.config \
  --genome GRCh38_hmf \
  --mode panel_resource_creation \
  --input samplesheet.panel_resource_creation.csv \
  --outdir output/ \
  \
  --driver_gene_panel DriverGenePanel.38.tsv \
  --target_regions_bed target_regions_definition.38.bed.gz \
  --isofox_gene_ids rna_gene_ids.csv # Optional, only provide if panel supports RNA data
```

### Overview of training procedure

The training procedure that `onconalyser` runs involves 2 mains steps which are described below.

#### Step 1: Run COBALT, AMBER, SAGE & ISOFOX

##### COBALT

COBALT is used determine raw read depth ratios (before GC normalisation).

**Commands**:

```
java -jar cobalt.jar \
    -tumor SAMPLE_ID \
    -tumor_bam /sample_data/SAMPLE_ID.bam \ 
    -output_dir /sample_data/cobalt \ 
    -gc_profile /ref_data/GC_profile.1000bp.37.cnp \
    -tumor_only_diploid_bed DiploidRegions.37.bed.gz \
    -threads 10 \ 
```
##### AMBER

AMBER is required to determine the gender of the samples for copy number normalisation. Note that we recommend to lower 
`-tumor_min_depth` to ensure sufficient coverage of heterozygous points.

**Commands**:

```
java -jar amber.jar com.hartwig.hmftools.amber.AmberApplication \
    -tumor SAMPLE_ID \
    -tumor_bam /sample_data/SAMPLE_ID.bam \
    -loci /path/to/GermlineHetPon.37.vcf.gz \ 
    -tumor_min_depth 2 \ 
    -output_dir /sample_data/ \
    -threads 10 \
```

##### SAGE

SAGE performs variant calling so that a panel-specific PON can be created.

**Commands**:

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

##### ISOFOX

Calculates gene and transcript expression with ISOFOX. ISOFOX is only run for panels supporting RNA sequencing data.

**Commands**:

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

> [!NOTE]
> ISOFOX requires the expected counts file to have been generated for the correct RNA read length (see [ISOFOX readme](https://github.com/hartwigmedical/hmftools/blob/master/isofox/README.md#transcript-expression)) for details.

#### Step 2: Run training normalisation scripts

##### COBALT: Target regions normalisation TSV

The target regions normalisation file contains the expected relative enrichment for each on target region. This file is created by COBALT 
with the following steps:
- For each 1K region covering any target region, extract each sample's tumor read count and the GC profile mappability and GC ratio bucket
- Calculate median and median read counts for each sample, and overall sample mean and median counts
- Normalise each sample's tumor read counts per region
- Calculate a median read count from all samples per GC ratio bucket
- Write a relative enrichment for each region to the output file, with a min enrichment of 0.1
- If no WGS is available for normalisation, the tumorGCRatio is assumed to be 1 for autosomes. The gender of each sample must be provided. 
Female samples are excluded from Y chromosome normalisation and males use a tumorGCRatio of 0.5 for the sex chromosomes

**Commands**:

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

##### PAVE: Panel artefact PON

In addition to the WGS PON, PAVE utilises a panel-specific PON to capture panel specific artefacts. Any non-hotspot variant found 3 or more 
times with a qual (TQP) of > 40 and modified map-qual factor of > -10 is added to the panel PON.

**Commands**:

```
java -cp pave.jar com.hartwig.hmftools.pave.pon_gen.PonBuilder \
  -sample_id_file training_sample_ids.csv \
  -vcf_path "/training_sample_data/*.sage.vcf.gz" \
  -ref_genome_version V38 \
  -output_dir /ref_data/ \
  -log_debug \
```

##### ISOFOX: TPM normalisation

This step is only required for panels with RNA coverage (e.g. TSO500) and is done to normalise the TPM so that it is equivalent to WTS. 
Specifically, The median adjusted TPM is calculated for each gene across the panel samples. If the median is zero, a replacement value of 
0.01 is used instead. The adjustment factor is calculated by dividing the panel median value by the corresponding whole genome value for 
each gene.

**Commands**:

```
java -cp isofox.jar com.hartwig.hmftools.isofox.cohort.CohortAnalyser \
  -root_data_dir /sample_data/ \
  -sample_data_file /sample_data/training_sample_ids.csv \
  -analyses PANEL_TPM_NORMALISATION \
  -gene_distribution_file /ref_data/isofox.wgs_gene_distribution.37.csv \
  -gene_id_file /ref_data/targeted_pane_gene_ids.csv \
  -output_dir /ref_data/ \
```

> [!NOTE] 
> The adjustment factors are calculated at the gene level and not at the transcript level. This means the adjusted TPMs for transcripts from 
> panel sequencing are not reliable.

## Targeted analysis pipeline differences 

### Functional differences per tool

<TO DO: make into a table>

#### COBALT
COBALT normalises copy number and masks off-target regions according to the CN normalisation file. If a targetRegions file is provided, then a target enrichment rate is calculated simply as the median tumorGCRatio for the specified regions.
Any depth windows outside of the targetRegions file are masked so that they are ignored downstream by PURPLE. Depth windows found in the TSV file are normalised first by the overall target enrichment rate for the sample, then by the relativeEnrichment for that depth window and finally by the normal GC bias adjustment. The GC bias is calculated using on target regions only.

#### AMBER
The following filters are applied:
* min_depth (in tumor) > 25
* Tumor ref and alt support >= 2
* Min_depth_percent and max_depth_percent are not applied
* Tumor ref and alt VAF >= 0.05
* AMBER loci must be within 300 bases of a target region

#### PURPLE
**MSI estimate**
To estimate MSI, a set of microsatellites with high coverage in the panel must also be defined.

For a set of microsatellite sites defined in the MSI target bed file count the number of passing variants at MSI sites ignoring SNV, MNV and 1 base deletes and requiring a VAF cutoff of > 0.15 for 2 and 3 base deletes or 0.08 for 4+ base deletes or any length insertion.

We estimate MSI rate as:
```
MSIndelsPerMb = 220 * # of MSI variants / # of MSI sites in panel
```

**TML & TMB estimate**
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


#### ISOFOX
TPM is normalised to bring panel gene expression in-line with WGS expression rates.

#### SAGE
SAGE is run in high depth mode. See SAGE readme for details.

### Recommended parameter values per tool
The following parameters are calibrated for panel sequencing and are set differently to WGS. These are the default panel parameter values.

#### AMBER

```
-target_regions_bed Target Regions BED file
```

#### COBALT
```
-target_region Target Regions Normalisation TSV
-pcf_gamma 50
```

#### SAGE
```
-high_depth_mode
```

#### PURPLE
```
-target_regions_bed Target Regions BED file
-target_regions_ratios Target Regions Ratios file
-target_regions_msi_indels Target Regions MSI Indels file
```

## Future improvements
* **HRD** - vCHORD (HRD classifier) is completed and will be made available in an upcoming release.
* **TOO** - vCUPPA (CNN based tissue of origin classifier) is under development
