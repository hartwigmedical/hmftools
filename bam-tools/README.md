# BamTools

Various routines to run on a BAM file

## BamMetrics
Capture BAM metrics

### Usage

```
java -cp bam-tools.jar com.hartwig.hmftools.bamtools.metrics.BamMetrics \
   -bam_file SAMPLE.bam \
   -ref_genome ref_genome.fasta \
   -output_dir /output_dir/ 
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
map_qual_threshold | Reads below this map quality count towards MAPQ counts, default 20
base_qual_threshold | Bases below this base quality count towards BASEQ count, default 10
max_coverage | Positions with coverage above this count towards CAPPED counts, default 250


## RegionSlicer
Slice a BAM and gather remote mate and supplementary reads

### Usage

```
java -cp bam-tools.jar com.hartwig.hmftools.bamtools.slice.RegionSlicer \
   -bam_file SAMPLE.bam \
   -ref_genome ref_genome.fasta \
   -ref_genome_version 37 \
   -output_dir /output_dir/ 
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
output_dir | Path for output file(s)
log_level | INFO or DEBUG
threads | Multi-thread count, default 1

## Version History and Download Links
- [1.0](https://github.com/hartwigmedical/hmftools/releases/tag/compar-v1.0)
