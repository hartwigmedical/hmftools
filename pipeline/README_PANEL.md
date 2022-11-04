# HMF DNA Targeted Example Pipeline

![HMF_Pipeline](./hmf_tools_panel_pipeline.png)

These scripts demonstrate how to run each HMF component in turn to produce DNA variant calling and analysis on a panel tumor BAM. 

They match the current tool version, configuration and resource files as used in the current HMF GCP pipeline (see [Platinum](https://github.com/hartwigmedical/platinum)). 

## Set-up

1. Download the latest release JAR for each tool as listed [here](https://github.com/hartwigmedical/hmftools#current-versions).
- also ensure that samtools (1.10 or higher) and bwa (0.7.17 or higher) are on the path

2. Download the resources files for either GRCh37 or GRCh38 from [HMFTools-Resources > DNA-Resources](https://console.cloud.google.com/storage/browser/hmf-public/HMFtools-Resources/dna_pipeline/).
- The latest resource files version is v5.31. 
- The latest resource files for the TSO-500 panel is labeled 'hmf_tso500_pipeline_resources.38_v5.31.gz'
- The reference genome files are available separately [HMFTools-Resources > Ref-Genome](https://console.cloud.google.com/storage/browser/hmf-public/HMFtools-Resources/ref_genome/).

3. Call the pipeline with the following arguments:
- a sample tumorId (eg 'COLO829T')
- the sample data directory with an existing directory named as per the sample's tumorId
- tumor BAM and BAM index files in the sample's directory, named as tumorId.bam
- all required tools in a tools directory
- all required resource files in a resource files directory
- the reference genome version - either 'V37' or 'V38'
- panel mode 'PANEL' (instead of 'WGS')
- number of threads used for each component
- maximum memory allocated to each component (default=12GB)

```
./scripts/run_pipeline ./scripts /sample_data/ /ref_data_dir/ /tools_dir/ COLO829T V38 PANEL 10 16 \
```  


