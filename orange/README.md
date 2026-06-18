# Oncogenic Results of Analyzing the Genome

Orange summarizes the key outputs from all algorithms in the Hartwig suite into a single PDF and JSON file.

Example report for a [whole genome sample](./docs/DemoWGS.orange.pdf).

Example report for a [targeted panel sample](./docs/DemoPanel.orange.pdf).

[Orange User Manual](./docs/orange_user_manual.pdf)


## Command

Orange uses most outputs from the WiGiTs pipeline.

### WGS mode

```
java -jar orange.jar \
    -experiment_type "PANEL"
    -tumor TUMOR_ID \
    -reference REFERENCE_ID \
    -primary_tumor_location Skin \
    -ref_genome_version V38 \
    -pipeline_sample_root_dir /sample_oa_results/
    -output_dir /path/to/where/to/write/output \
```

### Targeted Panel mode

```
java -jar orange.jar \
    -experiment_type "PANEL"
    -tumor TUMOR_ID \
    -panel_name "TSO500"
    -primary_tumor_location Skin \
    -pipeline_sample_root_dir /sample_oa_results/
    -output_dir /path/to/where/to/write/output \
```

### Arguments

| Argument                 | Description                                                                                                                   |
|--------------------------|-------------------------------------------------------------------------------------------------------------------------------|
| pipeline_version_file    | Path to the file containing the (platinum) pipeline version used.                                                             |
| sampling_date            | Optional. Sets the sampling date to the specified date. Expected format is YYMMDD. If omitted, sampling date will not be set. |
| analysis_date            | Optional. Sets the analysis date to the specified date. Expected format is YYMMDD. If omitted, current date is used as analysis date. |
| sequencing_type          | Illumina (default), SBX, Ultima                                                                                               |
| ref_genome_version       | V37 (default) or V38                                                                                                          |
| primary_tumor_location   | Printed at top of report                                                                                                      | 
| experiment_type          | WGS or PANEL                                                                                                                  |
| rna_sample_id            | Used to display RNA sample genotype values from SageAppend                                                                    |
| panel_name               | Optional, for display only                                                                                                    
| add_disclaimer           | If set, adds a "research use only" disclaimer to the footer of every page.                                                    |  
| pipeline_sample_root_dir | Optional, all individual algo paths are derived from this path, assuming the pipeline has been run using HMF pipeline         |
| sample_data_dir          | Optional, all data is expected to exist in the root of this path                                                              | 
