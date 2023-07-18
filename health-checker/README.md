# Health Checker

Checks Purple WC stats, Flagstats and BAM Metrics for key required values

## Usage

```
java -jar health-checker.jar \
   -tumor SAMPLE_T \
   -reference SAMPLE_R \
   -tum_flagstat_file SAMPLE_T.flagstats \
   -tum_wgs_metrics_file /SAMPLE_T.wgsmetrics \
   -ref_flagstat_file SAMPLE_R.flagstats \
   -ref_wgs_metrics_file /SAMPLE_R.wgsmetrics \
   -purple_dir /path_purple_files/ \
   -output_dir /output_dir/ 
```

## Arguments

Argument | Description 
---|---
tumor | Tumor sample ID
tum_flagstat_file / ref_flagstat_file | Sambamba or Samtools FlagStat file
tum_wgs_metrics_file / ref_wgs_metrics_file | BAM Metrics file as produced by Picard or HMF BamTools
purple_dir | Purple pipeline output files 
output_dir | Output directory
do_not_write_evaluation_file | If specified no output success/failure file is written


## Version History and Download Links
- [3.4](https://github.com/hartwigmedical/hmftools/releases/tag/health-checker-v3.4)
- [3.3](https://github.com/hartwigmedical/hmftools/releases/tag/health-checker-v3.3)
- [3.2](https://github.com/hartwigmedical/hmftools/releases/tag/health-checker-v3.2)
- [3.1](https://github.com/hartwigmedical/hmftools/releases/tag/health-checker-v3.1)
- [3.0](https://github.com/hartwigmedical/hmftools/releases/tag/health-checker-v3.0)
