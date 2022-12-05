# Oncogenic Results of Analyzing the Genome

ORANGE summarizes the key outputs from all algorithms in the Hartwig suite into a single PDF and JSON file:

1. The algo depends exclusively on config and data produced by the [Hartwig platinum pipeline](https://github.com/hartwigmedical/platinum)
   and hence can always be run as final step without any additional local data or config required.
2. The algo combines RNA and DNA data to present an integrated DNA/RNA analysis of a tumor sample.
   In general, ORANGE respects the mode in which the pipeline has been run (tumor-only, panel, whole genome etc).
3. ORANGE can be configured to convert all germline findings into somatic findings, thereby enabling downstream users to be oblivious
   of germline findings without missing anything actually relevant.
4. Everything that is labeled as a driver by any of the Hartwig algorithms is displayed in the PDF along with the driver likelihood.
   This effectively means that everything reported
   by [patient-reporter](https://github.com/hartwigmedical/oncoact/tree/master/patient-reporter)
   is present in the ORANGE pdf and json.
5. An additional exhaustive WGS and WTS scan is performed for anything interesting that may be potentially relevant but not picked up as a
   driver. Details of what is considered interesting are described in below.
6. A comprehensive range of QC measures and plots is displayed which provides in-depth details about the data quality of the tumor sample.

An example tumor + reference report based on the publicly available melanoma cell line COLO829 can be
found [here](src/main/resources/Test.orange.pdf).

Note that neither this readme nor the report itself contains any documentation about the Hartwig algorithms and output. For questions in
this area please refer to the specific algorithm documentation present
on [https://github.com/hartwigmedical/hmftools](https://github.com/hartwigmedical/hmftools)

The front page of the ORANGE report lists all high-level stats about the sample along with genome-wide visualisations of all mutations and
SNV/Indel clonality. In addition to this front page, the following chapters are present in the ORANGE report:

- [Somatic Findings](#somatic-findings): What potentially relevant mutations have been found in the tumor specifically?
- [Germline Findings](#germline-findings): What potentially relevant mutations have been found in the germline data?
- [Immunology](#immunology): What can we tell about the immunogenicity of the tumor sample?
- [RNA Findings](#rna-findings): What potentially relevant findings have we detected in RNA?
- [Cohort Comparison](#cohort-comparison): How do the various properties of this tumor compare to existing cancer cohorts?
- [Quality Control](#quality-control): Various stats and graphs regarding the quality of the data and interpretation thereof.

## Running ORANGE

### Base (tumor-only) mode

```
java -jar orange.jar \
   -tumor_sample_id tumor_sample
   -primary_tumor_doids doid1;doid2 \
   -ref_genome_version 37 \
   -output_dir /path/to/where/to/write/output \
   -doid_json /path/to/input_doid_tree.json \
   -cohort_mapping_tsv /path/to/input_cohort_mapping.tsv \
   -cohort_percentiles_tsv /path/to/input_cohort_percentiles.tsv \
   -driver_gene_panel_tsv /path/to/driver_gene_panel.tsv \
   -known_fusion_file /path/to/known_fusion_file.tsv \
   -tumor_sample_wgs_metrics_file /path/to/tumor_sample_wgs_metrics \
   -tumor_sample_flagstat_file /path/to/tumor_sample_flagstats \
   -sage_somatic_tumor_sample_bqr_plot /path/to/sage_tumor_sample_bqr_plot \
   -purple_data_directory /path/to/purple_data \
   -purple_plot_directory /path/to/purple_plots \
   -linx_somatic_data_directory /path/to/linx_somatic_data \
   -linx_plot_directory /path/to/linx_plots \
   -lilac_result_csv /path/to/lilac_results.csv \
   -lilac_qc_csv /path/to/lilac_qc.csv \
   -annotated_virus_tsv /path/to/annotated_virus.tsv \
   -chord_prediction_txt /path/to/chord_prediction.txt \
   -cuppa_result_csv /path/to/cuppa_results.tsv \
   -cuppa_summary_plot /path/to/cuppa_summary_plot \
```

### Additional parameters when running tumor-reference mode

```
    -reference_sample_id reference_sample \
    -ref_sample_wgs_metrics_file /path/to/reference_sample_wgs_metrics \
    -ref_sample_flagstat_file /path/to/reference_sample_flagstats \
    -safe_germline_gene_coverage_tsv /path/to/sage_germline_gene_coverage.tsv \
    -sage_somatic_ref_sample_bqr_plot /path/to/sage_ref_sample_bqr_plot \
    -linx_germline_data_directory /path/to/linx_germline_data \
    -peach_genotype_tsv /path/to/peach_genotypes.tsv 
```

### Additional parameters when RNA data is available

```
    -rna_sample_id rna_sample \
    -isofox_gene_distribution_csv /path/to/isofox_gene_distribution.csv \
    -isofox_alt_sj_cohort_csv /path/to/isofox_alt_sj_cohort.csv \
    -isofox_summary_csv /path/to/isofox_summary.csv \
    -isofox_gene_data_csv /path/to/isofox_gene_data.csv \
    -isofox_fusion_csv /path/to/isofox_fusion.csv \
    -isofox_alt_splice_junction_csv /path/to/isofox_alt_splice_junctions.csv
```

### Additional optional parameters across all modes

| Argument                    | Description                                                                                                                                        |
|-----------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------|
| pipeline_version_file       | Path to the file containing the (platinum) pipeline version used.                                                                                  |
| cuppa_feature_plot          | In case the cuppa summary does not fit onto one page, an additional cuppa feature plot is generated that has to be passed separately               |
| convert_germline_to_somatic | If set, converts all germline findings to somatic findings, thereby obfuscating the germline part of the analysis without loosing any actual data. |
| experiment_date             | Sets the experiment date to the specified date if set. Expected format is YYMMDD. If omitted, current date is used as experiment date.             |
| limit_json_output           | If set, limits all lists in the JSON output to a single entry to facilitate manual inspection of the JSON output.                                  |
| log_debug                   | If set, additional DEBUG logging is generated.                                                                                                     |

### Somatic Findings

In addition to all somatic drivers (SNVs/Indels, copy numbers, structural variants and fusions) the following is considered potentially
interesting and added to the report:

- Other potentially relevant variants:
    1. Variants that are hotspots or near hotspots but not part of the reporting gene panel.
    2. Exonic variants that are not reported but are phased with variants that are reported.
    3. Variants that are considered relevant for tumor type classification according to Cuppa.
    4. Variants with synonymous impact on the canonical transcript of a reporting gene but with a reportable worst impact
    5. Variants in splice regions that are not reported in genes with splice variant reporting enabled.
- Other regions with amps, or with deletions in other autosomal regions:
    1. Gains in genes for which we report amplifications with a relative minimum copy number between 2.5 and 3 times ploidy.
    2. Any chromosomal band location with at least one gene lost or fully amplified is considered potentially interesting.
       A maximum of 10 additional gains (sorted by minimum copy number) and 10 additional losses are reported as potentially interesting:
        - For a band with more than one gene amplified, genes are first picked based on heighest minimum copy number.
        - In case of a loss in a band with a reported loss, the additional loss is considered potentially interesting in case it is
          associated with clinical evidence.
        - For a band with a loss that has no losses reported in this band already, an arbitrary gene is picked.
- Other potentially relevant fusions. A maximum of 10 additional fusions (picked arbitrarily) are reported as potentially interesting:
    1. Any fusion that is not reported and has a reported type other than NONE.
    2. Any fusion in a gene that is configured as an oncogene in the driver gene panel.
- Other viral presence:
    1. Any viral presence that is not otherwise reported.
- Potentially interesting gene disruptions:
    1. Any unreported but disruptive gene disruption that is disrupting an exon which lies within a promiscuous exon range based on
       the fusion knowledgebase.
- Potentially interesting LOH events:
    1. In case MSI is detected, LOH (if present) is shown for the following genes: MLH1, MSH2, MSH6, PMS2, EPCAM
    2. In case HRD (based on CHORD) is detected, LOH (if present) is shown for the following genes: BRCA1, BRCA2, RAD51C, PALB2

In case ORANGE was run in DNA+RNA mode, DNA findings will be annotated with RNA:

- Driver and potentially interesting variants are annotated with RNA depth
- Driver and potentially interesting amps/dels are annotated with TPM, and corresponding percentile and foldChange for database and
  applicable tumor type
- Driver and potentially interesting fusions are annotated depending on fusion type:
    1. `EXON_DEL_DUP` and other intra-gene fusions are annotated with exon-skipping novel splice junctions
    2. @IG fusions are annotated with TPM of the 3' fusion gene
    3. Other fusions are annotated with RNA fusion details (detected fusions in RNA, and corresponding fragment support and depth of 5' and
       3' junction)

### Germline Findings

In addition to all germline SNV/Indel tumor drivers determined by [PURPLE](../purple), the following is added to the report:

- Other potentially relevant variants
    1. Any hotspots that are not configured to be reported.
    2. Any hotspots that are filtered based on quality.
- Missed variant likelihood (MVLH) per gene, presenting the likelihood of missing a pathogenic variant in case there would have been one
  present.
- Potentially pathogenic germline deletions
- Potentially pathogenic germline disruptions
- (Large-scale) germline CN aberrations.

Germline CN aberrations are determined by [PURPLE](../purple) and include aberrations such as klinefelter or trisomy X.

### Immunology

The immunology chapter is work-in-progress and will report on various immunology properties of the tumor sample.

The chapter currently presents the following:

- HLA-A/B/C details
    1. QC Status
    2. Detected alleles, annotated with #total fragments and somatic annotation (tumor copy number, #mutations)

In case ORANGE was run in DNA+RNA mode, the alleles will be annotated by RNA fragment support.

### RNA Findings

If run with RNA, this chapter displays potentially interesting RNA details:

- QC Details
- Drive gene panel genes with high TPM (>90th percentile database & tumor type) or low TPM (<5th percentile database or tumor type)
- Potentially interesting support for Known or Promiscuous fusions not detected in our DNA analysis pipeline
- Potentially interesting novel splice junctions
    1. Exon-skipping events in `EXON_DEL_DUP` fusion genes
    2. Novel exon/intron events in driver gene panel genes

### Cohort Comparison

The cohort comparison reports all the properties of a tumor sample that [Cuppa](../cuppa) considers for determining tumor type. The cohort
comparison displays the prevalence of the tumor's properties with respect to the cohorts that Cuppa could potentially assign the sample to:

- Genomic position distribution of SNVs and their tri-nucleotide signature
- Sample traits of the tumor (for example, number of LINE insertions)
- (Driver) features of the tumor.

Do note that RNA features and cohort comparison thereof are only included if ORANGE was run in combined DNA/RNA mode.

### Quality Control

The quality control chapter provides extensive details that can help with interpreting the overall [PURPLE](../purple) QC status or
investigate potential causes for QC failure.

- The high-level QC from [PURPLE](../purple)
- Various details from the tumor and reference samples flagstats and coverage stats
- Various plots from [PURPLE](../purple)
- BQR plots from both reference and tumor sample from [SAGE](../sage)

### Version History and Download Links

- [2.0](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v2.0)
    - Many and major datamodel changes in the ORANGE json output.
    - Remove PROTECT dependency from ORANGE including clinical evidence chapter and potentially interesting mutations based on evidence
    - Support tumor-only mode
    - Plots are now copied as part of ORANGE algo and paths to plots are registered relative to the location of the output json file.
    - Various minor changes and bug fixes to RNA:
        - Improve formatting of non-duplicate fragments in RNA findings
        - Germline variants are now actually annotated with RNA in case RNA data is present.
        - RNA sample config is mandatory in case RNA data is provided.
    - Disruptions are now displayed exactly as they are provided by linx without further conversions or merging.
- [1.10](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v1.10)
    - Classify each driver gene as wild-type in case:
        1. Purple fit is reliable
        2. No reportable events have been found in the gene.
    - Various additions and improvements to selection of potentially interesting events:
        - Unreported variants with canonical synonymous impact but reportable worst impact are added as potentially interesting variants.
        - Unreported variants in splice regions of genes for which we report splice variants are added as potentially interesting variants.
        - Unreported variants near hotspots are added as potentially interesting.
        - Unreported fusions of all oncogenes are added as potentially interesting fusions.
        - Driver genes that are nearly amplified (2.5 < minRelCN < 3) are added as potentially interesting gains.
        - Non-driver amplifications and losses are preferentially selected for being suspicious in case they have associated evidence.
        - Non-driver losses that have associated evidence but also a reported loss in the same band are retained as suspicious loss.
        - Breakends that are disruptive inside a pathogenic exon range of a promiscuous fusion gene but did not lead to fusions are added as
          potentially interesting disruptions.
    - Various new annotations of existing tables:
        - Gene disruptions are annotated by their cluster ID to be able to determine whether multiple events happened in same event.
            - New input: `linx_structural_variant_tsv`
        - Expression entries in the RNA chapter are annotated with their tumor CN
    - Clinical evidence is now grouped by event rather than treatment.
    - Various minor and technical improvements:
        - Add optional parameter `experiment_date` which, when provided in YYMMDD format, sets the experiment date. If param is not
          provided, the current date is used (as before).
        - Variants with identical phased inframe canonical effect are dedup'ed prior to reporting.
        - Treatment counts on the front page are retrieved from general evidence only and no longer include trials.
        - Average chromosome arm copy numbers are calculated and stored in ORANGE JSON:
            - New input: `purple_somatic_copy_number_tsv`
        - All potentially interesting sections are now also written to the JSON output.
        - The structural variant driver table has been removed as all interesting events where duplicated with events in other sections
        - Add optional parameter `limit_json_output`, which if set limits the size of the JSON output to facilitate manual inspection of
          JSON datamodel.
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
            - `isofox_summary_csv`, `isofox_gene_data_csv`, `isofox_fusion_csv`, `isofox_alt_splice_junction_csv` are the actual isofox data
              files.
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
    - Handle multiple drivers per gene where non-canonical transcripts are included. Current behaviour is to ignore non-canonical transcript
      drivers.
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
