# Qsee

Qsee collects and calculates QC metrics from pipeline outputs. Metrics not meeting warn/fail QC thresholds are flagged.


A visualisation (see [**example**](src/main/resources/COLO829v004T.qsee.vis.report.pdf)) is also generated providing an overview of all 
metrics and compares the respective sample to an existing cohort.

## Input files

Qsee accepts the below pipeline files as input:

- [**PURPLE**](https://github.com/hartwigmedical/hmftools/tree/master/purple)
  - `.purple.purity.tsv`: tumor only (mandatory)
  - `.purple.qc`: tumor only (mandatory)
- [**BamTools (BamMetrics)**](https://github.com/hartwigmedical/hmftools/tree/master/bam-tools#bammetrics)
  - `.bam_metric.summary.tsv`: tumor (mandatory), normal
  - `.bam_metric.coverage.tsv`: tumor (mandatory), normal
  - `.bam_metric.flag_counts.tsv`: tumor (mandatory), normal
  - `.bam_metric.frag_length.tsv`: tumor, normal
  - `.bam_metric.gene_coverage.tsv`: tumor, normal
- [**COBALT**](https://github.com/hartwigmedical/hmftools/tree/master/cobalt)
  - `.cobalt.gc.median.tsv`: tumor, normal
- [**ESVEE**](https://github.com/hartwigmedical/hmftools/tree/master/esvee)
  - `.esvee.prep.disc_stats.tsv`: tumor only
- [**REDUX**](https://github.com/hartwigmedical/hmftools/tree/master/redux)
  - `.redux.bqr.tsv`: tumor, normal
  - `.redux.duplicate_freq.tsv`: tumor only
  - `.redux.ms_table.tsv.gz`: tumor, normal

Only (some of) the PURPLE and BamTools files for the tumor sample are mandatory. This allows Qsee to run in tumor-only mode, as well as
to be flexible to missing input data.

## Running Qsee

### Single sample mode

1) Put all input files in a single directory (or provide [inputs in separate directories](#input-files-in-separate-directories)):

```
inputs/
├── TUMOR_1.cobalt.gc.median.tsv
├── REFERENCE_1.cobalt.gc.median.tsv
...
```

2) Run Qsee:

```
java -jar qsee.jar \
-tumor TUMOR_1 \
-reference REFERENCE_1 \
-sample_data_dir input/ \
-driver_gene_panel DriverGenePanel.38.tsv \
-cohort_percentiles_file qsee.cohort.percentiles.tsv.gz \
-output_dir output/
```

This will generate the following files:
- `TUMOR_ID.qsee.status.tsv.gz`: QC status (PASS/WARN/FAIL) per sample and per feature with a QC threshold defined
- `TUMOR_ID.qsee.vis.data.tsv.gz`: Visualisation raw data
- `TUMOR_ID.qsee.vis.report.pdf`: Visualisation comparing samples to an existing cohort, and highlighting samples not meeting QC thresholds (see [**example**](src/main/resources/COLO829v004T.qsee.vis.report.pdf))

Some key optional arguments:
- `-reference`: Can be omitted to run Qsee in tumor-only mode.
- `-cohort_percentiles_file`: If omitted, cohort percentiles will not be plotted. This is useful for panels with too few samples to generate cohort percentiles.
- `-targeted_mode`: Use default QC thresholds for targeted mode.
- `-threshold_overrides_file`: If provided, overrides the default thresholds for each feature (see [Threshold overrides](#threshold-overrides)).
- `-allow_missing_input`: If provided, Qsee will not fail if input files are missing. This is useful if for example only some BamMetrics 
  output files are available due to running an older pipeline version

### Multi-sample mode

Qsee can be run for multiple samples all at once.

1) Put all input files in a single directory (or provide [inputs in separate directories](#input-files-in-separate-directories)):

```
inputs/
├── TUMOR_1.cobalt.gc.median.tsv
├── REFERENCE_1.cobalt.gc.median.tsv
...
├── TUMOR_1.esvee.prep.disc_stats.tsv
├── REFERENCE_1.esvee.prep.disc_stats.tsv
...
├── TUMOR_1.purple.purity.tsv
├── TUMOR_1.purple.qc
...
├── TUMOR_1.bam_metric.summary.tsv
├── REFERENCE_1.bam_metric.summary.tsv
...
├── TUMOR_1.redux.ms_table.tsv.gz
├── REFERENCE_1.redux.ms_table.tsv.gz
...
```

2) Create a sample IDs file (in tumor-only mode, column ReferenceId can be omitted):

```
TumorId	 ReferenceId
TUMOR_1  REFERENCE_1
TUMOR_2  REFERENCE_2
TUMOR_3  REFERENCE_3
...      ...
```

3) Run Qsee:

```
java -jar qsee.jar \
-tumor TUMOR_ID \
-reference REFERENCE_ID \
-sample_data_dir pipeline_output/ \
-driver_gene_panel DriverGenePanel.38.tsv \
-cohort_percentiles_file qsee.cohort.percentiles.tsv.gz \
-output_dir output/ \
-threads 12 \
-merge_plots
```

This will generate the following files:
```
multisample.qsee.status.tsv.gz
multisample.qsee.vis.data.tsv.gz
multisample.qsee.vis.report.pdf
```

If `-merge_plots` is omitted, Qsee will generate a separate PDF for each sample.

### Threshold overrides

The threshold overrides file can be provided to arg `threshold_overrides_file` to override the default thresholds for each feature.

This file has the below format. `NaN` can be provided under column `Threshold` to set no threshold (i.e. always PASS).

```
SampleType  FeatureType    FeatureName        QcStatusType    ComparisonOperator  Threshold
TUMOR	    SUMMARY_TABLE  MAPPED_PROPORTION  FAIL            LESS_THAN           0.95
TUMOR       SUMMARY_TABLE  MEAN_COVERAGE      WARN            LESS_THAN           NaN
...
```

A default thresholds file with the above format can be written using the below command. This file shows all overridable thresholds
and can be used as a template for customising thresholds.

```
java -cp qsee.jar com.hartwig.hmftools.qsee.status.DefaultThresholdsWriter \
-output_dir output/ \
# -targeted_mode ## Add this arg to write default targeted mode thresholds
```

## Generating cohort percentiles

Qsee can compare samples to an existing cohort. The cohort percentiles file stores the percentiles for each metric for the tumor and/or
normal samples in the cohort.

As with [multi-sample mode](#multi-sample-mode), create a sample IDs file and put all input files in a single directory. Then run:

```
java -cp qsee.jar com.hartwig.hmftools.qsee.cohort.CohortPercentilesTrainer \
-sample_id_file sample_ids.txt \
-sample_data_dir inputs/ \
-driver_gene_panel DriverGenePanel.38.tsv \
-output_dir output/ \
-threads 12
```

Other optional arguments:
- `-write_cohort_features` If provided, write the raw feature values used to calculate the percentiles.
- `-allow_missing_input`: If provided, Qsee will not fail if input files are missing This is useful if e.g. only some BamMetrics files are
  available due to running an older pipeline version

## Input files in separate directories

For [Multi-sample mode](#multi-sample-mode) or [Generating cohort percentiles](#generating-cohort-percentiles), 
if input files are located in different tool directories, tool dir args can be used to provide those paths:

```
java -jar qsee.jar \
-sample_id_file sample_ids.txt \
-cobalt_dir cobalt/*/ \
-esvee_dir esvee/*/ \
-purple_dir purple/*/ \
-redux_tumor_dir redux/*/ \
-redux_ref_dir redux/*/ \
-bam_metrics_tumor_dir bamtools/*/ \
-bam_metrics_ref_dir bamtools/*/ \
... ## other required args
```

Wildcard (`*`) directories are optional and are substituted with the tumor IDs and normal IDs from the sample IDs file. 
Qsee will for example expect the REDUX MS table files at:
- `redux/TUMOR_1/TUMOR_1.redux.ms_table.tsv.gz`
- `redux/REFERENCE_1/REFERENCE_1.redux.ms_table.tsv.gz`

## All arguments

**Sample IDs**
- `-tumor`: Tumor ID
- `-reference`: Reference ID
- `-sample_id_file`: Sample ID TSV file containing columns TumorId and optionally ReferenceId

**Input/output**
- `-output_dir`: Output directory
- `-output_id`: Suffix to add to all output files
- `-sample_data_dir`: Directory containing all input files
- `-bam_metrics_tumor_dir`: Directory containing BamMetrics output for the tumor sample
- `-bam_metrics_ref_dir`: Directory containing BamMetrics output for the reference sample
- `-redux_tumor_dir`: Directory containing REDUX output for the tumor sample
- `-redux_ref_dir`: Directory containing REDUX output for the reference sample
- `-cobalt_dir`: Directory containing COBALT output
- `-esvee_dir`: Directory containing ESVEE output
- `-purple_dir`: Directory containing PURPLE output

**Resource files**
- `-driver_gene_panel`: Path to driver gene panel
- `-cohort_percentiles_file`: Path to cohort percentiles file (see [Generating cohort percentiles](#generating-cohort-percentiles))

**Other options**
- `-allow_missing_input`: If provided, Qsee will not fail if input files are missing
- `-merge_plots`: In multi-sample mode produce a single PDF with all samples.
- `-sequencing_type`: Can be `ILLUMINA`, `SBX` or `ULTIMA`. Only used to determine base qual bins in BaseQualRecalibrationPrep
- `-skip_categories`: A comma separated list of categories to skip prep for. Valid options: BASE_QUAL_RECALIBRATION, COVERAGE_DISTRIBUTION, DISCORDANT_FRAG_FREQ, DUPLICATE_FREQ, FRAG_LENGTH_DISTRIBUTION, GC_BIAS, MISSED_GENE_VARIANT, MS_INDEL_ERROR
- `-write_cohort_features`: If provided, write the raw feature values when calculating cohort percentiles. 
- `-threads`: Number of threads to use for parallel processing (one sample per thread)
- `-log_level`: Log level (default: `INFO`). See [log4j docs](https://logging.apache.org/log4j/2.x/manual/customloglevels.html) for valid values