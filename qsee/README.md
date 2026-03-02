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

The below example command (using arg `-sample_data_dir`) assumes all input files are in the same directory and runs Qsee in tumor/normal 
mode. For tumor-only mode, arg `-reference` can be omitted.

```
java -jar qsee.jar \
-tumor TUMOR_ID \
-reference REFERENCE_ID \
-sample_data_dir input/ \
-driver_gene_panel DriverGenePanel.38.tsv \
-cohort_percentiles_file qsee.cohort.percentiles.tsv.gz \
-output_dir output/
```

If input files are located in their respective tool directories, the below example command can be used:
```
java -jar qsee.jar \
-tumor TUMOR_ID \
-reference REFERENCE_ID \
-cobalt_dir cobalt/ \
-esvee_dir esvee/ \
-purple_dir purple/ \
-redux_tumor_dir redux/TUMOR_ID/ \
-redux_ref_dir redux/REFERENCE_ID/ \
-tumor_metrics_dir bamtools/TUMOR_ID/ \
-ref_metrics_dir bamtools/REFERENCE_ID/ \
-driver_gene_panel DriverGenePanel.38.tsv \
-cohort_percentiles_file qsee.cohort.percentiles.tsv.gz \
-output_dir output/
```

This will generate the following files:
- `TUMOR_ID.qsee.status.tsv.gz`: QC status (PASS/WARN/FAIL) per sample and per feature with a QC threshold.
- `TUMOR_ID.qsee.vis.report.pdf`: Visualisation comparing samples to an existing cohort, and highlighting samples not meeting QC thresholds.
- `TUMOR_ID.qsee.vis.data.tsv.gz`: Visualisation raw data

An example of the visualisation can be found here: [**COLO829v004T.qsee.vis.report.pdf**](src/main/resources/COLO829v004T.qsee.vis.report.pdf)

### Key optional arguments

Some arguments can be included or omitted for flexibility:
- `-reference`: Can be omitted to run Qsee in tumor-only mode.
- `-cohort_percentiles_file`: If omitted, cohort percentiles will not be plotted. This is useful for panels with too few samples to generate
  cohort percentiles.
- `-threshold_overrides_file`: If provided, overrides the default thresholds for each feature (see [Threshold overrides](#threshold-overrides)).
- `-allow_missing_input`: If provided, Qsee will not fail if input files are missing This is useful if e.g. only some BamMetrics files are 
  available due to running an older pipeline version

### Threshold overrides

The threshold overrides file can be provided to arg `threshold_overrides_file` to override the default thresholds for each feature. 

This file has the below format. `NaN` can be provided under column `Threshold` to set no threshold (i.e. always PASS).

```
SampleType  FeatureType    FeatureName        QcStatusType    ComparisonOperator  Threshold
TUMOR	    SUMMARY_TABLE  MAPPED_PROPORTION  FAIL            LESS_THAN           0.95
TUMOR       SUMMARY_TABLE  MEAN_COVERAGE      WARN            LESS_THAN           NaN
...
```

A default thresholds file can be written using the below command. This file shows all overridable thresholds, and can be used as a template
for customising thresholds.

```
java -cp qsee.jar com.hartwig.hmftools.qsee.status.DefaultThresholdsWriter -output_dir output/
```

## Generating cohort percentiles

Qsee can compare samples to an existing cohort. The cohort percentiles file stores the percentiles for each metric for the tumor and/or
normal samples in the cohort.

To generate this file, first create a sample IDs file (in tumor-only mode, column ReferenceId can be omitted).

```
TumorId	 ReferenceId
TUMOR_1  REFERENCE_1
TUMOR_2  REFERENCE_2
TUMOR_3  REFERENCE_3
...      ...
```

Then, given your pipeline outputs have the below structure:

```
pipeline_output/cobalt/
├── TUMOR_1.cobalt.gc.median.tsv
├── REFERENCE_1.cobalt.gc.median.tsv
...

pipeline_output/esvee/
├── TUMOR_1.esvee.prep.disc_stats.tsv
├── REFERENCE_1.esvee.prep.disc_stats.tsv
...

pipeline_output/purple/
├── TUMOR_1.purple.purity.tsv
├── TUMOR_1.purple.qc
...

pipeline_output/bamtools/
├──TUMOR_1.bam_metric.summary.tsv
├──REFERENCE_1.bam_metric.summary.tsv
...

pipeline_output/redux/
├──TUMOR_1.redux.ms_table.tsv.gz
├──REFERENCE_1.redux.ms_table.tsv.gz
...
```

the below command can be used:

```
java -cp qsee.jar com.hartwig.hmftools.qsee.cohort.CohortPercentilesTrainer \
-sample_id_file sample_ids.txt \
-cobalt_dir pipeline_output/cobalt/ \
-esvee_dir pipeline_output/esvee/ \
-purple_dir pipeline_output/purple/ \
-redux_tumor_dir pipeline_output/redux/ \
-redux_ref_dir pipeline_output/redux/ \
-tumor_metrics_dir pipeline_output/bamtools/ \
-ref_metrics_dir pipeline_output/bamtools/ \
-driver_gene_panel DriverGenePanel.38.tsv \
-output_dir output/ \
-threads 12
```

### Key optional arguments

- `-write_cohort_features` If provided, to write the feature values uses to calculate the percentiles.
- `-allow_missing_input`: If provided, Qsee will not fail if input files are missing This is useful if e.g. only some BamMetrics files are
  available due to running an older pipeline version

## All arguments

### Sample ID

| Argument             | Description                                                              |
|:---------------------|:-------------------------------------------------------------------------|
| `-tumor`             | Tumor ID                                                                 |
| `-reference`         | Reference ID                                                             |
| `-sample_id_file`    | Sample ID TSV file containing columns TumorId and optionally ReferenceId |

### Input/output

| Argument             | Description                                           |
|:---------------------|:------------------------------------------------------|
| `-output_dir`        | Output directory                                      |
| `-output_id`         | Suffix to add to all output files                     |
| `-sample_data_dir`   | Directory containing all input files                  |
| `-tumor_metrics_dir` | BamMetrics output dir for tumor sample                |
| `-ref_metrics_dir`   | BamMetrics output dir for reference sample            |
| `-redux_tumor_dir`   | REDUX output dir for tumor sample                     |
| `-redux_ref_dir`     | REDUX output dir for reference sample                 |
| `-cobalt_dir`        | COBALT output dir                                     |
| `-esvee_dir`         | ESVEE output dir                                      |
| `-purple_dir`        | PURPLE output dir                                     |
| `-sage_dir`          | SAGE somatic output dir (for backwards compatibility) |

### Resource files

| Argument                    | Description                                                                                                    |
|:----------------------------|:---------------------------------------------------------------------------------------------------------------|
| `-driver_gene_panel`        | Path to driver gene panel                                                                                      |
| `-cohort_percentiles_file`  | Path to [cohort percentiles file](#generating-cohort-percentiles) (not applicable to CohortPercentilesTrainer) |
| `-threshold_overrides_file` | Path to [threshold overrides file](#threshold-overrides) (not applicable to CohortPercentilesTrainer)          |

### Other options

| Argument               | Description                                                                                                                                                                                                                                              |
|:-----------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `-sequencing_type`     | Can be `ILLUMINA`, `SBX` or `ULTIMA`. Only used to determine base qual bins in BaseQualRecalibrationPrep                                                                                                                                                 |
| `-skip_categories`     | A comma separated list of categories to skip prep for. Valid options:<br>`BASE_QUAL_RECALIBRATION`, `COVERAGE_DISTRIBUTION`, `DISCORDANT_FRAG_FREQ`, `DUPLICATE_FREQ`,<br>`FRAG_LENGTH_DISTRIBUTION`, `GC_BIAS`, `MISSED_GENE_VARIANT`, `MS_INDEL_ERROR` |
| `-allow_missing_input` | If provided, Qsee will not fail if input files are missing                                                                                                                                                                                               |
| `-threads`             | Number of threads to use for parallel processing                                                                                                                                                                                                         |
| `-log_level`           | Log level (default: `INFO`). See [log4j docs](https://logging.apache.org/log4j/2.x/manual/customloglevels.html) for valid values                                                                                                                         |


 