# Crest - Check Reference Equality to Sample Transcriptome

To ensure that WTS samples are correctly matched to the same patient as
the WGS sample, Crest performs a simple test on a multi-sample VCF to ensure
that 90% of the SNPs have support in the specified RNA sample.

The input is assumed to be a germline VCF annotated with RNA calls. Thresholds
applied are described below and can be adjusted by the user.

The computed ratio of RNA supported to total reads is written to log output,
and a flag file "sample.CrestCheckSucceeded" or "sample.CrestCheckFailed"
is written to the output directory for use in multi-step pipelines.

All SNPS must have impact a gene with PASSED filters.

## Example usage

```bash
$ java -jar crest.jar -purple_dir /path/to/purple -sample COLO829v003T -rna_sample COLO829v003T_RNA
```

This assumes standard layout of the purple directory, with the wgs sample having been overwritten by
the sage annotated version. The purple/COLO829v003T.purple.germline.vcf.gz will be examined.

## Paramters

| Parameter                     | Description                                                                                     | Default               |
|-------------------------------|-------------------------------------------------------------------------------------------------|-----------------------|
| purple_dir                    | Location of annotated vcf                                                                       |                       |
| sample                        | Name of the WGS sample                                                                          |                       |
| rna_sample                    | If this path is set, all data is expected to exist in the root of this path                     |                       |
| do_not_write_evalutation_file | If given, the output .CrestCheck flag file is not produced                                      | false if not provided |
| min_total_reads               | Min number of reads at SNP in the RNA sample to count towards total                             | 10                    |
| min_rna_reads                 | Min number of reads at SNP matching the variant allele in RNA sample to count towards supported | 1                     |
| acceptance_ratio              | Lower threshold on ratio of rna supported / total reads for test to pass                        | 0.90                  |
| output_dir                    | Directory in which to write .CrestCheck flag file                                               |                       |

