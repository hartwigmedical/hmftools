# MarkDups

The HMF mark duplicates component performs both UMI aware and UMI agnositc duplicate marking. 
As the first component to run after alignment it also performs post-alignment improvements to the BAM, specifically by unmapping certain reads and deleting supplementary reads in specific problematic regions of the BAM

UMI are used to label each molecule in a sample with a unique sequence prior to PCR amplification.

The usage of UMIs is recommended primarily for three scenarios:  
* Counting of individual reads in low input samples
* very deep sequencing of RNA-seq libraries (> 80 million reads per sample),  
* detection of ultra-low frequency mutations in DNA sequencing.

UMI/Duplicate analysis is also a highly useful QC tool for library complexity and error rates

## Commands

```
java -jar mark-dups.jar 
    -sample SAMPLE_ID 
    -input_bam SAMPLE_ID.lane_01.bam,SAMPLE_ID.lane_02.bam,SAMPLE_ID.lane_03.bam  
    -ref_genome /path_to_fasta_files/
    -ref_genome_version V37
    -unmap_regions /ref_data/unmap_regions.37.tsv 
    -write_stats 
    -sambamba /path_to_sambamba/ 
    -samtools /path_to_samtools/ 
    -output_dir /path_to_output/
    -log_level DEBUG 
    -threads 24
    -multi_bam 
```

## Arguments

Argument | Required | Description
---|---|---
sample | Required | Sample ID
input_bam | Required | Path to BAM file(s)
output_bam | Optional | Output BAM file, otherwise will write SAMPLE_ID.mark_dups.bam
ref_genome | Required | Path to reference genome files as used in alignment
ref_genome_version | Required | V37 or V38
form_consensus | Optional | Form a consensus read from duplicates
unmap_regions | Optional | Regions of high depth, repeats or otherwise problematic for mapping
threads | Optional | Number of threads, default = 1
multi_bam | Optional | Write a BAM per thread prior to final merge, see note in 'Performance'
sambamba | Optional | Used to merge BAMs per thread when used with 'multi_bam' option and threads > 1
samtools | Optional | Used to sort and index final output BAM
output_dir | Optional | If not specified will write output same directory as input BAM
output_id | Optional | Additonal file suffix
read_output | Optional, default = NONE | Write detailed read info to CSV, types are: ALL, DUPLICATE, NONE
write_stats | Optional | Writes a duplicate frequency TSV file

### UMI Command

```
java -jar mark-dups.jar 
    -sample SAMPLE_ID 
    -input_bam SAMPLE_ID.aligned.bam
    -output_bam SAMPLE_ID.bam  
    -ref_genome /path_to_fasta_files/
    -ref_genome_version V37
    -unmap_regions /ref_data/unmap_regions.37.tsv 
    -umi_enabled
    -umi_duplex
    -umi_duplex_delim + 
    -umi_base_diff_stats 
    -sambamba /path_to_sambamba/ 
    -samtools /path_to_samtools/ 
    -output_dir /path_to_output/
    -log_level DEBUG 
    -threads 24
    -multi_bam 
```

### UMI Arguments

Argument | Description
---|---
umi_enabled | Extract UMI from Read ID and use in duplicate group identification
umi_duplex | Collapse duplex UMI groups
umi_duplex_delim | Duplex UMI character, default = '_'
umi_base_diff_stats | Write UMI statistics files


## Performance and Settings

When run wth multiple threads and 'multi_bam' enabled, a BAM will be written per thread and then merged and index at the end.
Recommended settings for a standard 100x tumor BAM is 16-24 CPUs and 48GB RAM.
Runtime on COLO829T with these settings is approximately 100mins.

## Output Files

File Name | Details 
---|---
SAMPLE_ID.duplicate_freq.tsv |  Frequency distribution of duplicate groups
SAMPLE_ID.reads.tsv | Detailed read information
SAMPLE_ID.umi_coord_freq.tsv | Frequency distribution of UMI duplicate groups
SAMPLE_ID.umi_edit_distance.tsv | Analysis of UMI differences and potential UMI base mismatches 
SAMPLE_ID.umi_nucleotide_freq.tsv | Frequency distribution for nucleotides in UMIs 

## Version History and Download Links
