# BamTools

Various routines to run on a BAM file

## BamMetrics
Capture BAM metrics

### Usage

```
java -cp bam-tools.jar com.hartwig.hmftools.bamtools.metrics.BamMetrics \
   -bam_file SAMPLE.bam \
   -ref_genome ref_genome.fasta \
   -output_dir /output_dir/ \
   -threads 10 \
```

### Configuration

Filter | Description
---|---
sample | Optional sample ID, used to set output filenames
bam_file | BAM file to analyse
ref_genome | Reference genome file used to create the BAM
ref_genome_version | Reference genome version, V37 (default) or V38
regions_bed_file | Optional BED file to restrict the sections of the genome analysed
specific_regions | Optional list of regions to analyse. Format 'chr:posStart-posEnd; etc'
write_old_style | Write a single output file with same format to Picard CollectWgsMetrics for compatibility
output_dir | Path for output file(s)
log_level | INFO or DEBUG
threads | Multi-thread count, default 1
partition_size | Default 1M bases, splits chromosomes to analyse in partitions
map_qual_threshold | Reads below this map quality count towards MAPQ counts, default 20
base_qual_threshold | Bases below this base quality count towards BASEQ count, default 10
max_coverage | Positions with coverage above this count towards CAPPED counts, default 250


## RegionSlicer
Slice a BAM and gather remote mate and supplementary reads

### Usage

```
java -cp bam-tools.jar com.hartwig.hmftools.bamtools.slice.RegionSlicer \
   -bam_file SAMPLE.bam \
   -output_prefix 
   -ref_genome ref_genome.fasta \
   -specific_regions 
   -output_dir /output_dir/ \
   -threads 10 \
```

### Configuration

Filter | Description
---|---
bam_file | BAM file to slice
output_prefix | Use in file output
ref_genome | Reference genome file used to create the BAM
regions_file | Optional BED file to restrict the sections of the genome analysed
partition_size | Default 1M bases, splits chromosomes to analyse in partitions
specific_regions | Optional list of regions to analyse. Format 'chr:posStart-posEnd; etc'
write_reads | Write reads to TSV file
drop_remote_supps | Ignore remote supplementary reads
output_dir | Path for output file, if omitted will write to same directory as BAM file
log_level | INFO or DEBUG
threads | Multi-thread count, default 1


## BamCompare
Compare 2 BAM/CRAM files and write differences into a TSV file.

### Usage

```
java -cp bam-tools.jar com.hartwig.hmftools.bamtools.compare.BamCompare \
   -orig_bam_file SAMPLE.original.bam \
   -new_bam_file SAMPLE.new.bam \
   -ref_genome ref_genome.fasta \
   -output_file SAMPLE.bam_compare.tsv.gz \ 
   -threads 10 \
```

### Configuration

| Filter                      | Required? | Description                                                                    |
|-----------------------------|-----------|--------------------------------------------------------------------------------|
| orig_bam_file               | Yes       | First BAM / CRAM file                                                          |
| new_bam_file                | Yes       | Second BAM / CRAM file                                                         |
| ref_genome                  | No        | Reference genome file, only required for CRAM input.                           |
| output_file                 | Yes       | Output comparison file (tsv or tsv.gz)                                         |
| ignore_dup_diffs            | No        | Ignore difference in duplicate flag                                            |
| ignore_alterations          | No        | Ignore consensus reads and internal unmappings                                 |
| ignore_consensus_reads      | No        | Ignore consensus reads                                                         |
| ignore_supplementary_reads  | No        | Ignore supplementary reads                                                     |
| partition_size              | No        | Default 10M bases, splits chromosomes to analyse in partitions                 |
| max_cached_reads_per_thread | No        | Maximum number of cached reads per thread, automatically calculated if not set |
| specific_regions            | No        | As above                                                                       |
| log_level                   | No        | WARN, INFO, DEBUG or TRACE                                                     |
| threads                     | No        | Multi-thread count, default 1                                                  |

## BamToFastq

Convert BAM / CRAM file into FASTQ

### Usage

```
java -cp bam-tools.jar com.hartwig.hmftools.bamtools.tofastq.BamToFastq \
   -bam_file SAMPLE.cram \
   -ref_genome ref_genome.fasta \
   -output_dir /output \
   -split_mode READ_GROUP \ 
   -threads 10 \
```

### Configuration

| Filter           | Required? | Description                                                                    |
|------------------|-----------|--------------------------------------------------------------------------------|
| bam_file         | Yes       | input BAM / CRAM file                                                          |
| ref_genome       | No        | Reference genome file, only required for CRAM input.                           |
| output_dir       | Yes       | Directory to write output into.                                                |
| output_id        | No        | Extra file name prefix for the output FASTQ files                              |
| split_mode       | No        | How output FASTQ files are split (NONE, READ_GROUP, THREAD)                    |
| specific_regions | No        | As above                                                                       |
| log_level        | No        | WARN, INFO, DEBUG or TRACE                                                     |
| threads          | No        | Multi-thread count, default 1                                                  |


## Version History and Download Links
- [1.0](https://github.com/hartwigmedical/hmftools/releases/tag/bam-tools-v1.0)
