# Oncogenic Results of Analyzing the Genome

ORANGE summarizes the key outputs from all algorithms in the Hartwig suite into a single PDF and JSON file (
see [orange-datamodel](../orange-datamodel)).

1. The algo depends exclusively on config and data produced by the [Hartwig genomic pipeline](../pipeline)
   and hence can always be run as final step without any additional local data or config required.
2. ORANGE respects the mode in which the pipeline has been run (tumor-only vs tumor-reference, targeted (panel) vs whole genome).
   In case RNA data is provided, the algo combines the RNA and DNA data to present an integrated DNA/RNA analysis of a tumor sample.
3. ORANGE can be configured to convert all germline driver variants to somatic driver variants, thereby obfuscating the germline driver part
   of the analysis without actually loosing this data.
4. Every event that is labeled as a driver by any of the Hartwig algorithms is displayed in the PDF along with the driver likelihood.
5. An additional exhaustive WGS and WTS scan is performed for anything interesting that may be potentially relevant but not picked up as a
   driver. Details of what is considered interesting are described in below.
6. A comprehensive range of QC measures and plots is displayed which provides in-depth details about the data quality of the samples
   provided.

Example reports based on the publicly available melanoma cell line COLO829 can be found here:

| Type                 | File                                                                                       | Note                                    | 
|----------------------|--------------------------------------------------------------------------------------------|-----------------------------------------|
| WGTS Tumor-reference | (please inquire at Hartwig)                                                                | We have no RNA data for COLO829         |
| WGS Tumor-reference  | [COLO829_WGS_TumorReference.pdf](src/main/resources/COLO829_WGS_TumorReference.orange.pdf) |                                         |
| Targeted Tumor-only  | [COLO829_Targeted_TumorOnly.pdf](src/main/resources/COLO829_Targeted_TumorOnly.orange.pdf) | (Derived from the tumor-reference data) |                  

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

Note that the JSON file contains every mutation found in the analysis and hence is much more extensive than the PDF.
The JSON file is meant to be used by downstream applications who wish to further interpret the results of the molecular analysis.

## Running ORANGE

ORANGE requires the output of various Hartwig algorithms, along with some resource files (doid json, cohort distributions, driver genes,
known fusions and ensembl data cache). The resource files required to run ORANGE can be
found [here](https://resources.hartwigmedicalfoundation.nl) for either 37 or 38 reference genome version.

### Base (targeted/panel) tumor-only DNA mode

```
java -jar orange.jar \
    -experiment_type "PANEL"
    -tumor_sample_id tumor_sample \
    -primary_tumor_doids "doid1;doid2" \
    -ref_genome_version "37" \
    -output_dir /path/to/where/to/write/output \
    -doid_json /path/to/input_doid_tree.json \
    -cohort_mapping_tsv /path/to/input_cohort_mapping.tsv \
    -cohort_percentiles_tsv /path/to/input_cohort_percentiles.tsv \
    -driver_gene_panel /path/to/driver_gene_panel.tsv \
    -known_fusion_file /path/to/known_fusion_file.tsv \
    -ensembl_data_dir /path/to/ensembl_data_directory \
    -tumor_sample_wgs_metrics_file /path/to/tumor_sample_wgs_metrics \
    -tumor_sample_flagstat_file /path/to/tumor_sample_flagstats \
    -sage_dir /path/to/sage_somatic_output \
    -purple_dir /path/to/purple_output \
    -purple_plot_dir /path/to/purple_plots \
    -linx_dir /path/to/linx_somatic_output \
    -linx_plot_dir /path/to/optional_linx_somatic_output_plots \
    -lilac_dir /path/to/lilac_output 
```

Note that `linx_plot_dir` is an optional parameter and can be left out completely in case linx has not generated any plots.

Note that `primary_tumor_doids` can be left blank (""). This parameter is used to look up cancer-type-specific percentiles for various
tumor characteristics. If primary tumor doids are not provided, percentiles are calculated against the full HMF database only.

### Additional parameters when whole genome tumor DNA data is available

```
   -virus_dir /path/to/virus_interpreter_output \
   -chord_dir /path/to/chord_output \
   -cuppa_dir /path/to/cuppa_output \
   -sigs_dir /path/to/sigs_output 
```

Also, the value of the `-experiment_type` parameter should be set to `WGS` for all whole genome configurations.

### Additional parameters when whole genome germline DNA data is available

```
    -reference_sample_id reference_sample \
    -ref_sample_wgs_metrics_file /path/to/reference_sample_wgs_metrics \
    -ref_sample_flagstat_file /path/to/reference_sample_flagstats \
    -sage_germline_dir /path/to/sage_germline_output \
    -linx_germline_dir /path/to/linx_germline_output \
    -peach_dir /path/to/peach_output.tsv 
```

### Additional parameters when whole genome RNA data is available

```
    -rna_sample_id rna_sample \
    -isofox_gene_distribution /path/to/isofox_gene_distribution.csv \
    -isofox_alt_sj_cohort /path/to/isofox_alt_sj_cohort.csv \
    -isofox_dir /path/to/isofox_output 
```

### Additional optional parameters across all modes

| Argument                    | Description                                                                                                                                                                                                                                                                           |
|-----------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| pipeline_version_file       | Path to the file containing the (platinum) pipeline version used.                                                                                                                                                                                                                     |
| sampling_date               | Sets the sampling date to the specified date if set. Expected format is YYMMDD. If omitted, current date is used as sampling date.                                                                                                                                                    |
| convert_germline_to_somatic | If set, converts all germline driver variants to somatic driver variants, thereby obfuscating the germline driver part of the analysis without actually loosing this data. Note that the data in other germline tables, except the pharmacogenetics table, is removed from this page. |
| add_disclaimer              | If set, adds a "research use only" disclaimer to the footer of every page.                                                                                                                                                                                                            |  
| limit_json_output           | If set, limits all lists in the JSON output to a single entry to facilitate manual inspection of the JSON output.                                                                                                                                                                     |
| log_debug                   | If set, additional DEBUG logging is generated.                                                                                                                                                                                                                                        |
| log_level                   | If set, overrides the default log level (INFO). Values can be `ERROR`, `WARN`, `INFO`, `DEBUG` and `TRACE`                                                                                                                                                                            |

### Additional run modes

Instead of individual algo directories, it is possible to configure a single directory in the following modes:

| Parameter                | Description                                                                                                                      | 
|--------------------------|----------------------------------------------------------------------------------------------------------------------------------|
| pipeline_sample_root_dir | If this path is set, all individual algo paths are derived from this path, assuming the pipeline has been run using HMF pipeline |
| sample_data_dir          | If this path is set, all data is expected to exist in the root of this path                                                      | 

## Somatic Findings

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
        - For a band with more than one gene amplified, the gene with the highest minimum copy number is picked.
        - For a band with a loss that has no losses reported in this band already, an arbitrary gene is picked.
- Other potentially relevant fusions. A maximum of 10 additional fusions (picked arbitrarily) are reported as potentially interesting:
    1. Any fusion that is not reported and has a reported type other than NONE.
    2. Any fusion in a gene that is configured as an oncogene in the driver gene panel.
- Other potentially interesting in-frame fusions in case no high drivers events are detected
    1. In case no high driver events are detected, any in-frame non chain terminated fusion that is not already reported  
- Other viral presence:
    1. Any viral presence that is not otherwise reported.
- Potentially interesting gene disruptions:
    1. Any unreported but disruptive gene disruption that is disrupting an exon which lies within a promiscuous exon range based on
       the fusion knowledgebase.
- Potentially interesting LOH events:
    1. In case MSI is detected, LOH (if present) is shown for the following genes: MLH1, MSH2, MSH6, PMS2, EPCAM
    2. In case HRD (based on CHORD) is detected, LOH (if present) is shown for the following genes: BRCA1, BRCA2, RAD51C, PALB2

In case ORANGE was run in DNA+RNA mode, DNA findings will be annotated with RNA:

- Drivers and potentially interesting variants are annotated with RNA depth
- Drivers and potentially interesting amps/dels are annotated with TPM, and corresponding percentile and fold change for database and
  applicable tumor type
- Drivers and potentially interesting fusions are annotated depending on fusion type:
    1. `EXON_DEL_DUP` and other intra-gene fusions are annotated with exon-skipping novel splice junctions
    2. @IG fusions are annotated with TPM of the 3' fusion gene
    3. Other fusions are annotated with RNA fusion details (detected fusions in RNA, and corresponding fragment support and depth of 5' and
       3' junction)

## Germline Findings

In addition to all germline SNV/Indel tumor drivers determined by [PURPLE](../purple), the following is added to the report:

- Other potentially relevant variants
    1. Any hotspots that are not configured to be reported.
    2. Any hotspots that are filtered based on quality.
- Potentially pathogenic germline deletions
- Potentially pathogenic germline LOH events
- Potentially pathogenic germline homozygous disruptions
- Potentially pathogenic germline gene disruptions
- Missed variant likelihood (MVLH) per gene, presenting the likelihood of missing a pathogenic variant in case there would have been one
  present.
- (Large-scale) germline CN aberrations.
    - Germline CN aberrations are determined by [PURPLE](../purple) and include aberrations such as klinefelter or trisomy X.
- Pharmacogenetics (DPYD & UGT1A1 status)

## Immunology

The immunology chapter reports on various immunology properties of the tumor sample.

The chapter presents the following:

- HLA-A/B/C details
    1. QC Status
    2. Detected alleles, annotated with #total fragments and somatic annotation (tumor copy number, #mutations)

In case ORANGE was run in DNA+RNA mode, the alleles will be annotated by RNA fragment support.

- Genetic immune escape analysis (inspired by [this paper](https://www.nature.com/articles/s41588-023-01367-1)). ORANGE attempts to detect the following mechanisms:
    - HLA-1 loss-of-function, detected in case one of the following mutations is present in either HLA-A, HLA-B or HLA-C:
        - MACN < 0.3 without the presence of a loss (proxy for LOH)
        - A clonal variant with canonical coding effect `NONSENSE_OR_FRAMESHIFT` or `SPLICE`
        - A clonal, biallelic variant with canonical coding effect `MISSENSE`
        - A full or partial loss
        - A homozygous disruption
    - Antigen presentation pathway inactivation, detected in case one of the following mutations is present in either B2M, CALR, TAP1, TAP2,
      TABBP, NLRC5, CIITA or RFX5:
        - A clonal variant with canonical coding effect `NONSENSE_OR_FRAMESHIFT` or `SPLICE`
        - A clonal, biallelic variant with canonical coding effect `MISSENSE`
        - A full or partial loss
        - A homozygous disruption
    - IFN gamma pathway inactivation, detected in case one of the following mutations is present in either JAK1, JAK2, IRF2, IFNGR1, IFNGR2,
      APLNR or STAT1
        - A clonal variant with canonical coding effect `NONSENSE_OR_FRAMESHIFT` or `SPLICE`
        - A clonal, biallelic variant with canonical coding effect `MISSENSE`
        - A full or partial loss
        - A homozygous disruption
    - (Potential) PD-L1 overexpression, detected in case CD274 is fully amplified.
    - CD58 inactivation, detected in case any of the following mutations happened in CD58:
        - A clonal variant with canonical coding effect `NONSENSE_OR_FRAMESHIFT` or `SPLICE`
        - A clonal, biallelic variant with canonical coding effect `MISSENSE`
        - A full or partial loss
        - A homozygous disruption
    - Epigenetics driven immune escape via SETDB1, detected in case SETDB1 is fully amplified.

## RNA Findings

If run with RNA, this chapter displays potentially interesting RNA details:

- QC Details
- Drive gene panel genes with high TPM (>90th percentile database & tumor type) or low TPM (<5th percentile database & tumor type)
- Potentially interesting support for known or promiscuous fusions not detected in our DNA analysis pipeline
- Potentially interesting novel splice junctions
    1. Exon-skipping events in `EXON_DEL_DUP` fusion genes
    2. Novel exon/intron events in driver gene panel genes

## Cohort Comparison

The cohort comparison reports all the properties of a tumor sample that [Cuppa](../cuppa) considers for determining tumor type. The cohort
comparison displays the prevalence of the tumor's properties with respect to the cohorts that Cuppa could potentially assign the sample to:

- Genomic position distribution of SNVs and their tri-nucleotide signature
- Sample traits of the tumor (for example, number of LINE insertions)
- (Driver) features of the tumor.

Do note that RNA features and cohort comparison thereof are only included if ORANGE was run in combined DNA/RNA mode.

## Quality Control

The quality control chapter provides extensive details that can help with interpreting the overall [PURPLE](../purple) QC status or
investigate potential causes for QC failure.

- The high-level QC from [PURPLE](../purple)
- Various details from the tumor and reference samples flagstats and coverage stats
- Various plots from [PURPLE](../purple)
- BQR plots from both reference and tumor sample from [SAGE](../sage)

## Version History and Download Links

- Upcoming:
  - Add unreported reason to fusions in ORANGE
  - Add etiology information to signatures and sort by allocation
  - Add percentage of unsupported segments to Quality control page
  - Show all viable fusions in ORANGE in samples where we detect no HIGH drivers
  - Ensure HIV is never reported in ORANGE report or included in ORANGE json
- [3.6.0](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v3.6.0)
    - Cuppa predictions with NaN likelihood are filtered in the ORANGE conversion of CUPPA results
- [3.5.1](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v3.5.1)
    - Workaround added for bug with mapping various ORANGE cohorts to non-existing ISOFOX cohorts.
- [3.5.0](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v3.5.0)
    - Add `PurpleTranscriptImpact.reported` and make `PurpleVariant.reported` a derived field. This
      uses the REPORTABLE_TRANSCRIPTS vcf field introduced in PURPLE 4.0.
    - Split pharmacogenetic haplotype into separate haplotype and genotype on front page and in table.
    - A new Genetic Immune Escape analysis has been added to datamodel and report. This analysis determines whether the sample uses any of
      the immune escape mechanisms known in cancer. See
      also [Genetic immune landscape paper](https://www.nature.com/articles/s41588-023-01367-1)
    - Technical: Removed work-around for pre v1.6 LILAC (incorrectly populating ref fragments in case of tumor-only)
- [3.4.0](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v3.4.0)
    - Produce ORANGE datamodel v2.4.0 (including max copy number for gene copy numbers)
    - Bugfix: ORANGE can now map stomach and esophageal squamous cell carcinomas to their rightful cohort.
        - Note: Available in ORANGE [3.3.1](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v3.3.1)
- [3.3.0](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v3.3.0)
    - Clonal likelihood is set to 0 or 1 based on variant copy number when converting germline variants to somatic
    - Bugfix: Fix bug in germline MVLH parsing that caused them to be underestimated by a factor 100.
    - Restrict germline MVLH table on PDF and related field in JSON to genes handled by SAGE germline.
- [3.2.0](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v3.2.0)
    - Support new ORANGE-datamodel (v2.3.0).
    - Support germline deletions heterozygous in the tumor:
        - Included in new section "Potentially pathogenic germline LOH events" on PDF report.
        - Included in new fields `allGermlineLossOfHeterozygosities` and `reportableGermlineLossOfHeterozygosities` in the JSON.
        - Converted to somatic LOH events when `convert_germline_to_somatic` parameter is provided.
        - Changed titles of somatic LOH and germline losses on PDF report to accommodate these changes.
    - Combine multiple germline loss calls for the same gene into one call.
    - Merge germline and somatic losses when both exist for the same gene.
    - Bugfix: Simple clusters affecting no exons are now excluded when counting expected number of linx plots
    - Bugfix: ORANGE throws an exception in case a cancer type is resolved for isofox that does not exist in the gene distribution data.
    - ORANGE throws an exception in case empty cuppa predictions are provided (cuppa output file is empty or is missing probabilities)
- [3.1.0](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v3.1.0)
    - ORANGE report now shows the new CUPPA v2.0 visualization
    - New ORANGE-datamodel (v2.2.0) supports the tabular output from both CUPPA v2.0 and CUPPA v1.x
    - `::` is used rather than `_` to concatenate the two genes fused in DNA and RNA fusions (impacting orange-datamodel!)
- [3.0.2](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v3.0.2)
    - Fixed bug with selecting high/low expressed genes for samples without a specific cancer-type (e.g. CUPs)
- [3.0.1](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v3.0.1)
    - Fixed bug (potential NPE) when resolving optional paths
- [3.0.0](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v3.0.0)
    - Disclaimer (if enabled) has bigger font size in footer, and a disclaimer is now also present in header.
    - New ORANGE-datamodel (v2.0.0) with lots of datamodel renames and clean-ups.
        - Parameter `experimentDate` has been renamed to `samplingDate`
        - Kataegis plots are now mandatory and expected to be produced in every situation by PURPLE.
    - Data is displayed in the report as "NA" in case of purple QC failure, in case the data by itself is not interpretable (e.g. TML).
    - The TMB status (high vs low) is displayed on the front page along with the actual TMB.
    - The status of UGT1A1 is displayed on the front page.
    - The undisrupted CN of a DUP, in case of a HOM_DUP_DISRUPTION, is displayed as undisrupted CN minus junction CN
    - Lilac RNA and ref counts are displayed as NONE in case they are not available (rather than 0)
    - The parameter `experiment_type` is now required, with valid values being PANEL or WGS.
    - Report all non-canonical variants in case of multiple variants in same gene
- [2.7.0](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v2.7.0)
    - Supports targeted tumor-only mode:
        - Omits Cuppa, Chord, Sigs and VirusBreakends
        - Omits Germline annotations
    - DOID cohort mapper throws exception instead of warn in case of invalid DOID combinations
    - Various DOID combinations are added to resolve to Esophagus or Stomach cohorts when combined with Gastroesophageal cancer DOIDs
    - A bug has been fixed with respect to using transcripts from the ensembl data cache that are not ensembl transcripts.
    - Test data and test report have been bumped to v5.33 pipeline
- [2.6.0](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v2.6.0)
    - Various updates to configuration:
        - All inputs are now configured via one directory per tool rather than individual files
        - `ensembl_data_directory` parameter has been renamed to `ensembl_data_dir`
        - `driver_gene_panel_tsv` parameter has been renamed to `driver_gene_panel`
        - `isofox_gene_distribution_csv` parameter has been renamed to `isofox_gene_distribution`
        - `isofox_alt_sj_cohort_csv` parameter has been renamed to `isofox_alt_sj_cohort`
        - `log_level` parameter has been added to allow manual override of the default log level
- [2.5.0](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v2.5.0)
    - Bugfix: Maintain linx clusters after ORANGE germline conversion
    - Include breakdown by classifier in CUPPA predictions
    - Fixed bug that used invalid RNA gene expression cohort percentiles for CUP tumors and prevented high-expression findings
    - Added new PURPLE variant effects: `PHASED_MISSENSE` and `PHASED_SYNONYMOUS`
- [2.4.1](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v2.4.1)
    - Fixes java.lang.IllegalStateException occurring when sample has virus integration.
- [2.4.0](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v2.4.0)
    - The ORANGE datamodel used in the json output has been separated from the ORANGE logic and available as an artifact for
      other projects to depend on (see also [ORANGE-datamodel](../orange-datamodel))
    - All copy numbers are rounded to single digit instead of no digits
    - Combination of urethra cancer and renal cell cancer is mapped to OTHER by ORANGE cohort mapper
    - The ORANGE cohort mapping application queries clinical view rather than datarequest
    - Added `-add_disclaimer` parameter that will print a "research use only" disclaimer in the footer when set
    - PDF documents now share a single instance of each font to reduce file size
    - Formatting for undetermined HRD type is improved
- [2.3](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v2.3)
    - Full support for reporting germline structural variants
        - This includes extended linx & purple datamodels with various germline findings
        - All somatic fields now have explicit "somatic" in their property name.
    - Added `-cuppa_chart_plot` parameter to hold the output of cuppa-chart
    - Cuppa, VirusInterpreter, Sigs and Chord are made optional to support pipeline running in targeted-mode
    - Add `experimentType` to ORANGE output which is either `TARGETED` or `FULL_GENOME` based on purple's targeted property.
    - Removed warning in case a variant potentially falls in the splice region of 2 neighbouring exons
- [2.2](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v2.2)
    - (COSMIC) signatures are added to the somatic findings (new parameter: `-sigs_allocation_tsv`)
    - Unreported germline variants are no longer converted to somatic variants in case "germline to somatic conversion" is enabled.
    - In case of NO_TUMOR, all germline variants are wiped in case "germline to somatic conversion" is enabled.
    - Populate affectedExon for intronic variants in splice regions.
    - RAD51B is added as a gene that is reported for LOH in case of HR deficiency
- [2.1](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v2.1)
    - Various improvements to RNA datamodel making the datamodel more explicit (enums rather than strings)
    - Virus interpreter data contains "all" and "reported" consistent with linx and purple.
    - affectedCodon and affectedExon are populated correctly in variant transcript impact (new parameter: `-ensembl_data_directory`)
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
