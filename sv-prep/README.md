# SV Prep - pre-GRIDSS filtering

SvPrep generates a filtered BAM file by identifying candidate SV reads, and then feeds this BAM into Gridss assembly.

## Running SvPrep

### Mandatory Arguments

Argument | Description 
---|---
sample | Sample ID
bam_file | Input BAM file
ref_genome | Reference genome fasta file
ref_genome_version | 37 (default) or 38
output_dir | Output directory
threads | Thread count

### Optional Arguments

Argument | Description 
---|---
known_fusion_bed | BED file with known fusion pair coordinates (as used in Gripss), require only 1 fragment for junctions
blacklist_bed | See below for explanation
existing_junction_file | Typically used for reference sample after tumor has been run - ensure fragment support is captured for these junctions
write_types | From list of 'JUNCTIONS', 'READS', 'BAM' and 'FRAGMENT_LENGTH_DIST', default is JUNCTIONS and BAM
partition_size | Default is 1 million bases
min_align_bases | Min required aligned bases for a junction fragment, default = 50
min_junction_frags | Min fragments to call a junction, default = 2

```
java -jar sv-prep.jar 
  -sample SAMPLE_ID
  -bam_file /sample_data/SAMPLE_ID.bam
  -ref_genome /path_to_ref_genome_fasta/
  -ref_genome_version [37 or 38] 
  -ensembl_data_dir /path_to_ensembl_files/
  -output_dir /sample_data/output/ 
```


## Overview and algorithm

### Junction Fragments


### Support Fragments


### Discordang Group Fragments


### Blacklist Regions



# Version History and Download Links
- [1.0](https://github.com/hartwigmedical/hmftools/releases/tag/sv-prep-v1.0)
