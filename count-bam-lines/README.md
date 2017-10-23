# COBALT 1.2

**Co**unt **ba**m **l**ines is designed to count the number of read starts within each 1000 base window of a bam.

It will only include reads that match ALL the following criteria:
* Equals or exceeds min quality (default 10)
* Is not unmapped
* Is not duplicated
* Is neither secondary nor supplementary

It then applies GC Normalization to calculate the read ratios. It can optionally apply a diploid normalization (to be used for reference samples)

Finally, the Bioconductor copy number package is used to generate segments from the ratio file.


## R Dependencies
Segmentation is done with the Bioconductor [copynumber](http://bioconductor.org/packages/release/bioc/html/copynumber.html) package.

This can be installed in R with the following commands:
```
   source("https://bioconductor.org/biocLite.R")
   biocLite("copy number")
```


## Usage

Argument | Default | Description
---|---|---
-reference_bam | None | Location of reference bam file. Can be omitted to regenerate ratios if prior cobalt file already exists.
-tumor_bam | None | Location of tumor bam file. Can be omitted to regenerate ratios if prior cobalt file already exists.
-output_dir | None | Directory to write output
-reference | None | Name of reference sample
-tumor | None | Name of tumor sample
-threads | 4 | Number of threads to use
-window_size | 1000 | Size of window
-min_quality | 10 | Min quality
-gc_profile | None | Location of GC profile.

Arguments without default values are mandatory.

### Example Usage

```
java -jar cobalt.jar -reference REFERENCE -reference_bam /run_dir/REFERENCE.bam -tumor TUMOR -tumor_bam /run_dir/TUMOR.bam -output_dir /run_dir/cobalt -threads 24 -gc_profile /path/to/GC_profile.1000bp.cnp
```


## Output
The following tab delimited files are written:

`/run_dir/cobalt/TUMOR.cobalt`

`/run_dir/cobalt/TUMOR.chr.len`

`/run_dir/cobalt/TUMOR.cobalt.ratio.pcf`

`/run_dir/cobalt/REFERENCE.cobalt.ratio.pcf`

TUMOR.cobalt contains the counts and ratios of the reference and tumor.

TUMOR.cobalt.ratio.pcf contains the segmented regions determined from the ratios.