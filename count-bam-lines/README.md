# COBALT 1.1

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
-input | None | Location of bam file. Can be omitted to regenerate ratios if prior count file already exists.
-output_dir | None | Directory to write output
-sample | None | Name of sample
-diploid | N/A | Applies diploid normalization. Used for reference samples.
-threads | 4 | Number of threads to use
-window_size | 1000 | Size of window
-min_quality | 10 | Min quality
-gc_profile | None | Location of GC profile.

Arguments without default values are mandatory.

### Example Usage

```
java -jar cobalt.jar -input /run_dir/SAMPLE.bam -output_dir /run_dir/cobalt -sample SAMPLE -threads 24 -gc_profile /path/to/GC_profile.1000bp.cnp
```


## Output
The following tab delimited files are written:

`/run_dir/cobalt/SAMPLE.cobalt.count`

`/run_dir/cobalt/SAMPLE.chr.len`

`/run_dir/cobalt/SAMPLE.cobalt.ratio`

`/run_dir/cobalt/SAMPLE.cobalt.ratio.pcf`

SAMPLE.cobalt.count contains the Chromosome, Position and ReadCount of that window.
Positions start from 1. A count of -1 indicates no read starts within that window.
Not all empty windows will be written to file - only the first and last of each chromosome.

SAMPLE.chr.len is written with the lengths of each chromosome from the bam file.
This assists with downstream processes to produce accurate window sizes.

SAMPLE.cobalt.ratio contains the Chromosome, Position and Ratio of each window.

SAMPLE.cobalt.ratio.pcf contains the segmented regions determined from the ratios.