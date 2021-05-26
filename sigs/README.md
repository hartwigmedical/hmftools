# Sigs
Fits sample SNV counts to trinucleotide signature definitions.

## Configuration

To install, download the latest compiled jar file from the [download links](#version-history-and-download-links). 

### Resource Files
Sigs requires a signature definition file, for example the COSMIC 30 SNV signatures. 
These can be downloaded here [HMFTools-Resources](https://resources.hartwigmedicalfoundation.nl/):

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
position_bucket_size | 0 | Default is to not calculate position-bucket frequencies
max_sample_count | 20000 | Cap any individual position bucket's SNV count to this level 

### Example Usage

```
java -jar sigs.jar \
   -sample COLO829T 
   -signatures_file /reference_file/snv_cosmic_signatures.csv \
   -somatic_vcf_file /sample_data/COLO829T.purple.somatic.vcf.gz \
   -output_dir /output_dir/ \
```
