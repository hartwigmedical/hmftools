# Crest - Check Reference Equality to Sample Transcriptome

To ensure that WTS samples are correctly matched to the same patient as
the WGS sample, Crest performs a simple test on a multi-sample VCF to ensure
that 90% of the germline SNPs have support in the specified RNA sample.

The input is assumed to be a germline VCF annotated with RNA calls. The
thresholds applied are described below and can be adjusted by the user.
Only SNPs impacting a gene with filters PASSED are counted.

The computed ratio of RNA supported to total reads is written to log output,
and a flag file "{sample}.CrestCheckSucceeded" or "{sample}.CrestCheckFailed"
is written to the output directory for use in multi-step pipelines.

## Example usage

```bash
$ java -jar crest.jar -purple_dir /path/to/purple -sample COLO829v003T -rna_sample COLO829v003T_RNA
```

This assumes standard layout of the purple directory, with the wgs sample having been overwritten by
the sage annotated version. The vcf file purple/COLO829v003T.purple.germline.vcf.gz is assumed
to exist and will be examined.

## Parameters

| Parameter         | Description                                                                                | Default               |
|-------------------|--------------------------------------------------------------------------------------------|-----------------------|
| purple_dir        | Location of annotated vcf                                                                  |                       |
| sample            | Name of the WGS sample, used to construct the VCF filename                                 |                       |
| rna_sample        | The name of the RNA sample in the vcf to be examined                                       |                       |
| do_not_write_file | If given, the output .CrestCheck flag file is not produced                                 | false if not provided |
| min_total_reads   | Min number of reads at SNP in the RNA sample to count towards total                        | 10                    |
| min_rna_reads     | Min number of reads at SNP matching the variant allele in RNA sample to count as supported | 1                     |
| acceptance_ratio  | Lower threshold on ratio of rna supported / total reads for test to pass                   | 0.90                  |
| output_dir        | Directory in which to write .CrestCheck flag file                                          |                       |

