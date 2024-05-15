# Sigs
Fits sample SNV counts to trinucleotide signature definitions.

## Configuration

To install, download the latest compiled jar file from the [download links](#version-history-and-download-links). 

### Resource Files
Sigs requires a signature definition file, for example the COSMIC 30 SNV signatures. 
These can be downloaded here  [HMFTools-Resources > DNA Pipeline](https://console.cloud.google.com/storage/browser/hmf-public/HMFtools-Resources/dna_pipeline/):

### Mandatory Arguments

Argument | Description 
---|---
sample | Sample ID
signatures_file | SNV signature definitions file
somatic_vcf_file | SNV VCF file
output_dir | Output directory for sample SNV fits

### Optional Arguments

Argument | Default | Description 
---|---|---
min_alloc| 0.01 | Allocate counts to a signature if exceeds this number of total SNV counts
min_alloc_perc| 0.0005 | Allocate counts to a signature if exceeds this % of total SNV counts
 
### Position Frequency Arguments
Sigs can generate genomic position frequencies by rounding positions to a specified position bucket size. 

Argument | Default | Description 
---|---|---
position_bucket_size | 0 | 0 means no position frequencies are generated, otherwise bucket size
max_sample_count | 20000 | Cap any individual position bucket's SNV count to this level 

### Example Usage

VCF sourced:
```
java -jar sigs.jar \
   -sample COLO829T 
   -signatures_file /reference_file/snv_cosmic_signatures.csv \
   -somatic_vcf_file /sample_data/COLO829T.purple.somatic.vcf.gz \
   -output_dir /output_dir/ \
```

# Version History and Download Links
- [1.2](https://github.com/hartwigmedical/hmftools/releases/tag/sigs-v1.2)
- [1.1](https://github.com/hartwigmedical/hmftools/releases/tag/sigs-v1.1)
- [1.0](https://github.com/hartwigmedical/hmftools/releases/tag/sigs-v1.0)
