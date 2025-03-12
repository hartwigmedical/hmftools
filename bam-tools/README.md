# BamTools

Various routines to run on a BAM file:
* [Bam Metrics](#bammetrics)
* [Region Slicer](#regionslicer)
* [Bam to Fastq](#bamtofastq)
* [Bam Compare](#bamcompare)
* [Bam Checker](#bamchecker)


## BamMetrics
Captures and writes the follow metrics for a BAM file:
- overall statistics for types of reads and coverage
- fragment length distribution
- coverage distribution
- flag stats - frequency by read flag type
- partition read stats per 1M (default) bucket
- targeted panel coverage if configured

### Usage

```
java -cp bam-tools.jar com.hartwig.hmftools.bamtools.metrics.BamMetrics \
   -bam_file SAMPLE.bam \
   -ref_genome ref_genome.fasta \
   -output_dir /output_dir/ \
   -threads 16 \
```

### Configuration

Filter | Description
---|---
sample | Optional sample ID, used to set output filenames
bam_file | BAM file to analyse
ref_genome | Reference genome file used to create the BAM
ref_genome_version | Reference genome version, V37 (default) or V38
output_dir | Path for output file(s)
log_level | INFO (default) or DEBUG
threads | Multi-thread count, default 1

### Optional Configuration

Filter | Description
---|---
regions_file | Optional BED file to restrict the sections of the genome analysed
only_target | Only capture metrics within the defined regions
specific_regions | Optional list of regions to analyse. Format 'chr:posStart-posEnd; etc'
partition_size | Default 1M bases, splits chromosomes to analyse in partitions
map_qual_threshold | Reads below this map quality count towards MAPQ counts, default 20
base_qual_threshold | Bases below this base quality count towards BASEQ count, default 10
max_coverage | Positions with coverage above this count towards CAPPED counts, default 250
exclude_zero_coverage | Exclude regions of zero coverage from overall statistics


## RegionSlicer
Slice a BAM across 1 or more regions, and gather any remote mate and supplementary reads linked to those regions

### Usage

```
java -cp bam-tools.jar com.hartwig.hmftools.bamtools.slice.RegionSlicer \
   -bam_file SAMPLE.bam \
   -output_prefix 
   -ref_genome ref_genome.fasta \
   -specific_regions chr10:89000000-91000000 \
   -bamtool /tools/sambamba \
   -output_dir /output_dir/ \
   -threads 16 \
```

### Configuration

Specify either a set of specific regions or pass in a regions TSV or BED file: 

Filter | Description
---|---
bam_file | BAM file to slice
output_prefix | Use in file output
ref_genome | Reference genome file used to create the BAM
regions_file | BED or TSV file to restrict the sections of the genome analysed
specific_regions | List of regions to analyse. Format 'chr:posStart-posEnd; etc'
output_dir | Path for output file, if omitted will write to same directory as BAM file
log_level | INFO (default) or DEBUG
threads | Multi-thread count, default 1

### Optional Configuration

Filter | Description
---|---
partition_size | Default 1M bases, splits chromosomes to analyse in partitions
bamtool | Tool for sorting and indexing output BAM (samtools or sambamba)
write_reads | Write reads to TSV file
drop_excluded | Ignores remote reads in the poly-G region on chr 2


## BamCompare
Compare 2 BAM/CRAM files and write differences into a TSV file.

### Usage

```
java -cp bam-tools.jar com.hartwig.hmftools.bamtools.compare.BamCompare \
   -orig_bam_file SAMPLE.original.bam \
   -new_bam_file SAMPLE.new.bam \
   -ref_genome ref_genome.fasta \
   -output_file SAMPLE.bam_compare.tsv.gz \ 
   -threads 16 \
```

### Configuration

Filter  | Description                                         
---|----
orig_bam_file      | First BAM / CRAM file
new_bam_file       | Second BAM / CRAM file
ref_genome         | Reference genome file, only required for CRAM input 
output_file        | Output comparison file (tsv or tsv.gz)
log_level         | As above
threads           | As above

### Optional Configuration

Filter | Description                                                    
---|---
ignore_dup_diffs   | Ignore difference in duplicate flag                            
ignore_alterations | Ignore consensus reads and internal unmappings                 
ignore_consensus_reads | Ignore consensus reads                                         
ignore_supp_reads  | Ignore supplementary reads
ignore_supp_attribute | Ignore differences in supplementary attribute values           
partition_size     | Default 10M bases, splits chromosomes to analyse in partitions
max_cached_reads_per_thread | Maximum number of cached reads per thread, automatically calculated if not set
specific_regions   | As above


## BamToFastq

Convert BAM / CRAM file into FASTQ. Drops any consensus reads added by Redux.

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

 Filter     | Description
------------|---
bam_file   | input BAM / CRAM file
ref_genome | Reference genome file, only required for CRAM input
output_dir | Directory to write output into
output_id  | Extra file name prefix for the output FASTQ files
log_level | As above
threads    | As above

### Optional Configuration

Filter | Description                                                           
---|----
split_mode       | How output FASTQ files are split (NONE, READ_GROUP (default), THREAD
specific_regions | As above


## BamChecker
Checks a BAM for fragment consistency and sets Mate CIGAR if missing from paired reads

### Usage

```
java -cp bam-tools.jar com.hartwig.hmftools.bamtools.checker.BamChecker \
   -bam_file SAMPLE.bam \
   -ref_genome ref_genome.fasta \
   -bamtool /tools/sambamba \
   -write_incompletes \
   -output_dir /output_dir/ \
   -threads 16 \
```

### Configuration

Filter | Description
---|---
bam_file | BAM file to slice
ref_genome | Reference genome file used to create the BAM
output_dir | Path for output file, if omitted will write to same directory as BAM file
log_level | As above
threads | As above

### Optional Configuration

Filter | Description
---|---
specific_regions | List of regions to analyse. Format 'chr:posStart-posEnd; etc'
bamtool | As above
write_incompletes | Write information about incomplete or invalid fragments to TSV file



## AltContigRemapper

Remaps HLA alt contig alignments using a supplied reference genome.

Supplementary or secondary HLA alignments are skipped.

### Usage

```
java -cp bam-tools.jar com.hartwig.hmftools.bamtools.remapper.AltContigRemapper \
   -orig_bam_file sample_with_hla_alt_contigs.bam \
   -ref_genome ref_genome.fasta \
   -output_file remapped.bam \
   -bamtool /opt/tools/samtools/1.20/samtools \
   -threads 6
```


### Configuration

| Filter        | Required? | Description                                                          |
|---------------|-----------|----------------------------------------------------------------------|
| orig_bam_file | Yes       | input BAM / CRAM file                                                |
| ref_genome    | Yes       | Reference genome file, used to remap alt contigs.                    |
| output_file   | Yes       | File into which the records (remapped and unaltered) are written.    |
| bamtool       | No        | Path to samtools executable. If supplied, the output will be sorted. |
| threads       | No        | Multi-thread count, default 1                                        |

