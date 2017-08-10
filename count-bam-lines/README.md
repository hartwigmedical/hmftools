# COBALT

**Co**unt **ba**m **l**ines is designed to count the number of read starts within each 1000 base window of a bam.

## Usage

Argument | Default | Description
---|---|---
-input | None | Location of bam file
-output_dir | None | Directory to write output
-sample | None | Name of sample
-threads | 4 | Number of threads to use
-window_size | 1000 | Size of window

Arguments without default values are mandatory.

### Example Usage

```
java -jar cobalt.jar -input /run_dir/SAMPLE.BAM -output_dir /run_dir/cobalt -sample SAMPLE -threads 24
```

## Output
The output is a tab delimited file with the name of the sample followed by the extension "cobalt.count".
Each row contains the Chromosome, Position and ReadCount of that windown. Positions start from 1.

A count of -1 indicates no read starts within that window.
Not all empty windows will be written to file - only the first and last of each chromosome.

A second output file "sample.chr.len" is written with the lengths of each chromosome from the bam file.
This assists with downstream processes to produce accurante last window sizes.