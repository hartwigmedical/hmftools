# Oncogenic Results of Analyzing the Genome

ORANGE summarizes the key outputs from all algorithms in the Hartwig suite into a single PDF and JSON file:
 1. The algo depends exclusively on config and data produced by the [Hartwig platinum pipeline](https://github.com/hartwigmedical/platinum) 
 and hence can always be run as final step without any additional local data or config required. 
 1. The algo combines RNA and DNA data to present an integrated DNA/RNA analysis of a tumor sample. 
 1. Everything that is labeled as a driver by any of the Hartwig algo's is displayed in the PDF along with the driver likelihood. 
 This effectively means that everything reported by [patient-reporter](../patient-reporter) is present in the ORANGE pdf and json.
 1. An additional exhaustive WGS and WTS scan is performed for anything interesting that may be potentially relevant but not picked up as a driver.
 Details of what is considered interesting are described in below.
 1. A comprehensive range of QC measures and plots is displayed which provides in-depth details about the data quality of the tumor sample. 
 
An example report based on the publicly available melanoma cell line COLO829 can be found [here](src/main/resources/Test.orange.pdf).

Note that neither this readme nor the report itself contains any documentation about the Hartwig algorithms and output. For questions in 
this area please refer to the specific algorithm documentation present on [https://github.com/hartwigmedical/hmftools](https://github.com/hartwigmedical/hmftools)

The front page of the ORANGE report lists all high-level stats about the sample along with genome-wide visualisations of all mutations and 
SNV/Indel clonality. In addition to this front page, the following chapters are present in the ORANGE report:
 
  - [Somatic Findings](#somatic-findings): What potentially relevant mutations have been found in the tumor specifically?
  - [Germline Findings](#germline-findings): What potentially relevant mutations have been found in the germline data? 
  - [Immunology](#immunology): What can we tell about the immunogenicity of the tumor sample?
  - [RNA Findings](#rna-findings): What potentially relevant findings do we detect in RNA?
  - [Cohort Comparison](#cohort-comparison): How do the various properties of this tumor compare to existing cancer cohorts?
  - [Clinical Evidence](#clinical-evidence): What genomic evidence has been found in favor of, or against, specific treatments?
  - [Quality Control](#quality-control): Various stats and graphs regarding the quality of the data and interpretation thereof. 

### Optional Configuration

Argument | Description
---|---
disable_germline  | If set, disables the germline findings chapter and transforms germline variants to somatic variants.
max_evidence_level | If set, filters evidence down to this level. For example, if "B" is passed as a parameter, only treatments with at least A or B level evidence are displayed in the clinical evidence chapter of the report. Do note that the front page always lists the count of all evidence present, regardless of this filter setting.

### Somatic Findings 

In addition to all somatic drivers (SNVs/Indels, copy numbers, structural variants and fusions) the following is considered potentially
interesting and added to the report:
 - Other potentially relevant variants
    1. Variants that are hotspots but not part of the reporting gene panel.
    1. Variants which have clinical evidence but are not part of the reporting gene panel.
    1. Coding variants that are not reported but are phased with variants that are reported.
    1. Variants that are considered relevant for tumor type classification according to Cuppa.
 - Other regions with amps, or with deletions in other autosomal regions:
    1. Any chromosomal band location with at least one gene lost or fully amplified is considered potentially interesting.
        - For a band with more than one gene amplified, the gene with the highest minimum copy number is picked.
        - For a band with a loss that has no losses reported in this band already, a random gene is picked.
    1. A maximum of 10 additional gains (sorted by minimum copy number) and 10 additional losses are reported as potentially interesting. 
 - Other potentially relevant fusions:
    1. Any fusion that is not reported and has a reported type other than NONE is picked. 
    1. Any fusion with clinical evidence is picked. 
    1. A maximum of 10 additional fusions (randomly picked) are reported as potentially interesting.
 - Other viral presence
    1. Any viral presence that is not otherwise reported is reported as potentially interesting.
 - Potentially interesting LOH events
    1. In case MSI is detected, LOH (if present) is shown for the following genes: MLH1, MSH2, MSH6, PMS2, EPCAM
    1. In case HRD (based on CHORD) is detected, LOH (if present) is shown for the following genes: BRCA1, BRCA2, RAD51C, PALB2

In case ORANGE was run in DNA+RNA mode, DNA findings will be annotated with RNA:
 - Driver and potentially interesting variants are annotated with RNA depth
 - Driver and potentially interesting amps/dels are annotated with TPM, and corresponding percentile and foldChange for database and applicable tumor type
 - Driver and potentially interesting fusions are annotated depending on fusion type:
    1. `EXON_DEL_DUP` and other intra-gene fusions are annotated with exon-skipping novel splice junctions
    1. @IG fusions are annotated with TPM of the 3' fusion gene
    1. Other fusions are annotated with RNA fusion details (detected fusions in RNA, and corresponding fragment support and depth of 5' and 3' junction)  
    
### Germline Findings

In addition to all germline SNV/Indel tumor drivers determined by [PURPLE](../purple), the following is added to the report:
 - Other potentially relevant variants
    1. Any hotspots that are not configured to be reported.
    1. Any hotspots that are filtered based on quality. 
 - Missed variant likelihood (MVLH) per gene, presenting the likelihood of missing a pathogenic variant in case there would have been one present.
 - Potentially pathogenic germline deletions
 - Potentially pathogenic germline disruptions
 - (Large-scale) germline CN aberrations.

Germline CN aberrations are determined by [PURPLE](../purple) and include aberrations such as klinefelter or trisomy X. 

### Immunology

The immunology chapter is work-in-progress and will report on various immunology properties of the tumor sample. 

The chapter currently presents the following:
- HLA-A/B/C details
    1. QC Status
    1. Detected alleles, annotated with #total fragments and somatic annotation (tumor copy number, #mutations)
    
In case ORANGE was run in DNA+RNA mode, the alleles will be annotated by RNA fragment support. 

### RNA Findings

If run with RNA, this chapter displays potentially interesting RNA details:
 -  QC Details
 -  Drive gene panel genes with high TPM (>90th percentile database & tumor type) or low TPM (<5th percentile database or tumor type)
 -  Potentially interesting support for Known or Promiscuous fusions not detected in our DNA analysis pipeline
 -  Potentially interesting novel splice junctions
    1. Exon-skipping events in `EXON_DEL_DUP` fusion genes
    1. Novel exon/intron events in driver gene panel genes

### Cohort Comparison

The cohort comparison reports all the properties of a tumor sample that [Cuppa](../cuppa) considers for determining tumor type. The cohort
comparison displays the prevalence of the tumor's properties with respect to the cohorts that Cuppa could potentially assign the sample to:
 - Genomic position distribution of SNVs and their tri-nucleotide signature
 - Sample traits of the tumor (for example, number of LINE insertions)
 - (Driver) features of the tumor.
 
Do note that RNA features and cohort comparison thereof are only included if ORANGE was run in combined DNA/RNA mode.  
 
### Clinical Evidence 
 
 The following algo is used to render clinical evidence in the ORANGE report based on [PROTECT](../protect) output:
  1. Evidence is split up based on applicable and "potentially interesting" based on PROTECT reported yes/no.
  1. Evidence is split between trials and non-trials which are further split up based on on/off label. 
  1. Evidence is grouped by treatment and split up between responsive and resistance evidence.
  1. Evidence is filtered based on the optional `max_reporting_level` configuration. 

### Quality Control

The quality control chapter provides extensive details that can help with interpreting the overall [PURPLE](../purple) QC status or 
investigate potential causes for QC failure.
 - The high-level QC from [PURPLE](../purple)
 - Various details from the tumor and reference samples flagstats and coverage stats
 - Various plots from [PURPLE](../purple)  
 - BQR plots from both reference and tumor sample from [SAGE](../sage)

### Version History and Download Links
- [1.9](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v1.9)
  - Proper support for RNA
    - Add (mandatory) `driver_gene_panel_tsv` and `known_fusion_file` inputs to support interpretation of isofox results
    - RNA Depth for variants is picked up in case purple somatic/germline variants have been annotated with RNA. 
    - Amps and dels are annotated with expression data (tpm, plus percentiles and fold change)
    - Fusions are annotated by RNA support:
        - For `EXON_DEL_DUP` and other intra-gene fusions, a list of novel splice junctions is used for annotation.
        - IG fusions are annotated with the expression data of the 3' gene
        - All other fusions are annotated with their equivalent counterparts in RNA.
    - A new chapter is added with RNA statistics and various types of novel findings that were not found in DNA.
    - Config impact:
        - `rna_sample_id` enables RNA-annotated variant loading when configured as-expected.
        - `isofox_gene_distribution_csv` and `isofox_alt_sj_cohort_csv` enable annotation of isofox results.
        - `isofox_summary_csv`, `isofox_gene_data_csv`, `isofox_fusion_csv`, `isofox_alt_splice_junction_csv` are the actual isofox data files. 
  - (Technical) improvements to CUPPA
    - CUPPA data loader favors overall combined score > DNA combined score > RNA combined score.
    - CUPPA data loader retains the combined prediction for every cancer type, not just the best prediction.
  - Add DPYD status on front page
  - `ref_genome_version` configures the ref genome version used, and is propagated into JSON output and report (QC chapter).
  - Support for germline SVs
    - `linx_germline_disruption_tsv` configures the path towards the LINX germline disruptions.
    - Reported germline disruptions are displayed in the Germline Findings chapter.
  - Support for germline deletions
    - `purple_germline_deletion_tsv` configures the path towards the PURPLE germline deletions.
    - Reported germline deletions are displayed in the Germline Findings chapter.
  - Support for LILAC
    - `lilac_result_csv` and `lilac_qc_csv` configure the paths to the LILAC data files
    - Immunology chapter displays the LILAC QC status along with the HLA alleles found.
  - Addition of potentially relevant LOH events:
    - In case of HRD: LOH is reported for BRCA1, BRCA2, PALB2, RAD51C
    - In case of MSI: LOH is reported for MLH1, MSH2, MSH6, PMS2, EPCAM
  - Evidence labeled as having `NO_BENEFIT` from PROTECT are filtered out for reporting.
- [1.8](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v1.8)
  - Only show every source name once in clinical evidence section.
  - Fix bug with selecting variants that are unreported but have evidence.
  - Improve support for multiple drivers on same gene but with different transcripts.
- [1.7](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v1.7)
  - Support for PROTECT v2.1
  - Support multiple LPS per variant in SAGE
- [1.6](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v1.6)
  - Transform germline variants to somatic in case germline is switched off (somatic findings + drivers on front page)
  - Add "upstream" to variant details in case variant is upstream without annotation
  - Add driver likelihood for viruses
  - Generify trial sources to include any trial source that is labeled as trial by SERVE
- [1.5](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v1.5)
  - Support for [PAVE](../pave/README.md)
  - Handle multiple drivers per gene where non-canonical transcripts are included. Current behaviour is to ignore non-canonical transcript drivers.
  - Support for proper fusion table rendering in case fusions are wrapped over multiple pages. 
- [1.4](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v1.4)
  - Fix a formatting problem in clinical evidence in case of very long genomic events
  - Support display of variants with coding impact relative to the 3' UTR region of a gene. 
  - Fix another bug with percentile determination in the absence of a cancer type
- [1.3](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v1.3)
  - Fix bug for generating reports for samples without a known cancer type
  - Fix bug when passing a Cuppa feature plot image that does not exist. 
- [1.2](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v1.2) 
  - Support for [Virus Interpreter v1.1](../virus-interpreter/README.md), including addition of % covered, mean coverage and expected 
  clonal coverage
  - The Cuppa best prediction is always displayed on the front page regardless of reliability of prediction.
  - More details about HR deficiency are displayed on front page in case sample is HR deficient
  - Pan-cancer and cancer-type specific percentiles for SV TMB are displayed on the front page
  - Other autosomal regions with deletions are no longer filtered for germline events, so all autosomal deletions are now displayed
  regardless of whether they occurred in germline or not.
  - Many technical and param changes described in linked release notes 
- [1.1](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v1.1)
  - Add JSON output of comprehensive platinum output
- [1.0](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v1.0)
  - Initial release
