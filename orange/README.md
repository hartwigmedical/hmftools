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
5. A comprehensive range of QC measures and plots is displayed which provides in-depth details about the data quality of the samples
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
    -tumor tumor_sample \
    -primary_tumor_doids "doid1;doid2" \
    -ref_genome_version "37" \
    -doid_json /path/to/input_doid_tree.json \
    -driver_gene_panel /path/to/driver_gene_panel.tsv \
    -pipeline_sample_root_dir /sample_oa_results/
    -tumor_metrics_dir /sample_oa_results/bamtools/tumor/
    -ref_metrics_dir /sample_oa_results/bamtools/reference/
    -tumor_redux_dir /sample_oa_results/alignments/tumor/
    -ref_redux_dir /sample_oa_results/alignments/reference/
    -output_dir /path/to/where/to/write/output \
```

Note that `linx_plot_dir` is an optional parameter and can be left out completely in case linx has not generated any plots.

Note that `primary_tumor_doids` can be left blank (""). This parameter is used to look up cancer-type-specific percentiles for various
tumor characteristics. If primary tumor doids are not provided, percentiles are calculated against the full HMF database only.

Also, the value of the `-experiment_type` parameter should be set to `WGS` for all whole genome configurations.

### Additional parameters when whole genome RNA data is available

```
    -rna_sample_id rna_sample \
    -isofox_dir /path/to/isofox_output 
```

### Additional optional parameters across all modes

| Argument                    | Description                                                                                                                                                                                                                                                                           |
|-----------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| pipeline_version_file       | Path to the file containing the (platinum) pipeline version used.                                                                                                                                                                                                                     |
| sampling_date               | Sets the sampling date to the specified date if set. Expected format is YYMMDD. If omitted, current date is used as sampling date.                                                                                                                                                    |
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

## Immunology

The immunology chapter reports on various immunology properties of the tumor sample.

The chapter presents the following:

- HLA-A/B/C details
    1. QC Status
    2. Detected alleles, annotated with #total fragments and somatic annotation (tumor copy number, #mutations)

In case ORANGE was run in DNA+RNA mode, the alleles will be annotated by RNA fragment support.

- Genetic immune escape analysis (inspired by [this paper](https://www.nature.com/articles/s41588-023-01367-1)). ORANGE attempts to detect
  the following mechanisms:
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
- Support for known or promiscuous fusions not detected in our DNA analysis pipeline
- Reportable novel splice junctions

## Cohort Comparison

The cohort comparison reports all the properties of a tumor sample that [Cuppa](../cuppa) considers for determining tumor type. The cohort
comparison displays the prevalence of the tumor's properties with respect to the cohorts that Cuppa could potentially assign the sample to:

- Genomic position distribution of SNVs and their tri-nucleotide signature
- Sample traits of the tumor (for example, number of LINE insertions)
- (Driver) features of the tumor.

Do note that RNA features and cohort comparison thereof are only included if ORANGE was run in combined DNA/RNA mode.

## Quality Control

See [QSee](../qsee) 
