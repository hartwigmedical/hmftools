# Redux

The Redux component performs both UMI aware and UMI agnostic duplicate marking. 
As the first component to run after alignment it also performs post-alignment improvements to the BAM, specifically by unmapping certain reads and deleting supplementary reads in specific problematic regions of the BAM

UMI are used to label each molecule in a sample with a unique sequence prior to PCR amplification.

The usage of UMIs is recommended primarily for three scenarios:  
* Counting of individual reads in low input samples
* very deep sequencing of RNA-seq libraries (> 80 million reads per sample),  
* detection of ultra-low frequency mutations in DNA sequencing.

UMI/Duplicate analysis is also a highly useful QC tool for library complexity and error rates

## Commands

```
java -jar redux.jar 
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
```

## Arguments

Argument | Required | Description
---|---|---
sample | Required | Sample ID
input_bam | Required | Path to BAM file(s)
output_bam | Optional | Output BAM file, otherwise will write SAMPLE_ID.redux.bam
ref_genome | Required | Path to reference genome files as used in alignment
ref_genome_version | Required | V37 or V38
form_consensus | Optional | Form a consensus read from duplicates
unmap_regions | Optional | Regions of high depth, repeats or otherwise problematic for mapping
threads | Optional | Number of threads, default = 1
sambamba | Optional | Used to merge BAMs per thread when used with threads > 1
samtools | Optional | Used to sort and index final output BAM
output_dir | Optional | If not specified will write output same directory as input BAM
output_id | Optional | Additonal file suffix
read_output | Optional, default = NONE | Write detailed read info to CSV, types are: ALL, DUPLICATE, NONE
write_stats | Optional | Writes a duplicate frequency TSV file

### UMI Command

```
java -jar redux.jar 
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
```

### UMI Arguments

Argument | Description
---|---
umi_enabled | Extract UMI from Read ID and use in duplicate group identification
umi_duplex | Collapse duplex UMI groups
umi_duplex_delim | Duplex UMI character, default = '_'
umi_base_diff_stats | Write UMI statistics files

## Algorithm

### 'Unmapping’ in problematic regions 

Due to the incompleteness and repetitiveness of the human genome and idiosyncratic errors in sequencing, there are certain regions that have reads aligned to them, often with very high coverage, which recurrently cause issues in downstream tools. Key examples include PolyG regions, long PolyA/T repeats, long dinucleotide repeats, ribosomal regions on short arms of Chr21 and chr22 and many centromeric regions.  1 key case is that duplicate fragments with one read consisting mainly of homopolymers will not be properly marked as duplicates as the homopolymer read may be mapped to numerous places in the genome 

We therefore introduce logic to handle problematic regions of recurrent high depth across germline samples and long repeats (see appendix for definition of problematic regions). Specifically, we unmap a read if it has <10 bases aligned outside a problematic region AND meet at least one of these criteria: 
- Region has very high depth (3rd highest maxDepth >1000, ~700kb in hg38) 
- Read pair is discordant – an inversion, translocation, one read unmapped or fragment length > 1000 bases 
- Read has at least 20 base of soft clip 

All supplementary and secondary reads with <10 bases aligned outside a problematic region are also deleted. 

Note that when a read is unmapped or a supplementary is deleted, other reads in the read group pair are also updated to reflect the mates unmmaped status.  Removing / unmapping these reads simplifies and improve performance of variant calling downstream including in SAGE, COBALT and SV calling. 

### Deduplication

There are 2 steps in the deduplication algorithm:

#### 1. Identify duplicates 
Duplicates can be marked either with or without UMIs. If UMIs are not available, then simply mark any fragments (including any secondary & supplementary alignments) with the same fragment strand orientation and identical unclipped 5’ alignment coordinates (allowing for strand) on both reads in fragment are identified as belonging to the same duplicate group.  For non-paired reads the length of the fragment must also be identical. For paired reads, the MC (mate cigar) tag is used, where available (note BWA normally populates this but not STAR), to identify from the 1st read whether the fragment is in the same duplicate group. Where it is not available, both reads are assessed together, at the cost of higher runtime memory. If a primary read is marked as duplicate, all supplementary and secondary reads from the same fragment will also be marked duplicate. If one of the reads is unmapped then all fragments with the same start coordinates on the mappable side are placed in the same duplicate group. 

For each duplicate group, the fragment with the highest R1 average base qual (or highest base qual  are assessed together) is kept as primary and the remaining are marked as duplicates. If UMIs are available, then 1 duplicate group is made per UMI allowing for a small edit distance using a directional network (as described by [UMI-tools](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5340976/)).  To implement this, UMIs with the same alignment coordinates are ranked from most supported to least supported.   Starting with the 2nd highest supported distinct UMI, each UMI may be collapsed recursively into a more supported UMI group if it is within 1 edit distance of any UMI in that group.   Following this, a 2nd single level sweep is made from most supported to least supported, where any UMIs with same alignment coordinates and a configurable edit distance [3].  A 3rd single level sweep is to a default edit distance of 4 to collapse any highly imbalanced groups (ie where 1 group has more than 25x the support of another). 

For DUPLEX UMIs only, duplicates with the same coordinates but opposite orientations are further collapsed into a single group and and marked as 'DUAL_STRAND’. 

If run in duplex mode, the UMI pairs must match (within an edit distance of 1) between the fragments to be collapsed. If multiple UMI groups with the same coordinates may be collapsed then the group with the highest support on the forward strand is matched first to the group with the highest support on the reverse strand  

#### 2. Consensus fragments 

Optionally, a new ‘consensus’ read/fragment can be created to represent each duplicate group. In this case all fragments in a duplicate group will be marked as a duplicate and the new consensus read will be the primary alignment. 

To construct the consensus fragment, the following logic is applied separately for read 1 and read 2: 
- The consensus cigar is chosen as the cigar with the most support followed by the least soft clip and then arbitarily. 
- Using the consensus cigar as the reference,  for each base assess each reads in the duplicate group and set the consensus base to the most supported base by sum of base qual. If 2 or more alleles have the same base qual support choose the reference first. 
- For reads which differ by an indel from the consensus, ignore the differences, and ensure that the subsequent bases match the position in the consensus cigar.  Ignore any soft clip base that is not present in the consensus cigar. 
- Set the base qual = max(supportedBaseQual) * [max(0,Sum(supportedBaseQual) - Sum(contraryBaseQual))] / Sum(supportedBaseQual) 

In the case of DUPLEX UMIs, the logic is applied to each strand individually and then again to merge the 2 strands.   When merging strands, if the consensus base on each strand is different and one matches the ref genome then the base is set to the ref genome, otherwise the highest qual base is chosen as per above.   

The ‘CR’ flag is added to the bam to record the number of reads contributing to each consensus read. The ‘UT’ tag is also used to mark the UMI group as either ‘SINGLE’ or ‘DUAL_STRAND’ 

## Performance and Settings

When run wth multiple threads, a BAM will be written per thread and then merged and index at the end.
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

### Duplicate Frequency

Field | Description 
---|---
DuplicateReadCount  | # of duplicates in group 
Frequency | Total read observations 
DualStrandFrequency  | Total dual strand observations 

### UMI Edit Distance 
This file shows the edit distance from each duplicate ead to each consensus UMI for coordinates with 1 or 2 distinct UMI groups only and can be used to detect potential underclustering of UMIS 

Field | Description 
---|---
UMICountWithCoordinatesField | # of unique UMI groups that share the coordinates (either 1 or 2) 
DuplicateReadCount | Total # of duplicates at coordinates 
Frequency  | Total read observations 
<ED_0 to ED_10>  | Count of reads by edit distance 

### UMI coordinate frequency

This file shows how frequently distinct UMI groups that share the same start position also share the same coordinates.    

Field | Description 
---|---
UniqueCoordinatesWithStartPosition | # of fragments with unique coordinates that share the start positoin 
PrimaryFragmentsWithStartPosition  | # of fragments after deduplication that share the start position 
Frequency  | Total read observations 

### UMI nucleotide frequency  

This file gives the nucleotide frequency by UMI base position 

Field | Description 
---|---
UMIPosition | Base position in UMI 
ACount | Frequency of ‘A’ NT 
CCount  | Frequency of ‘C’ NT 
GCount  | Frequency of ‘G’ NT 
TCount  | Frequency of ‘T’ NT 
NCount | Frequency of ‘N’ NT 

## Problematic regions file defintion

Certain regions of the genome are consistently problematic and lead to frequently obvious mismapped regions which can cause several downstream intepretation problems in variant calling.   We specifically identify 2 types of such regions - very long repeats and regions which consistently align with much higher coverage across many WGS samples and create a bed file   

Problematic high depth regions are identified (for both 37 and 38) by bserving the depth of 10 WGS germline samples with depth of approximate 40x.    
- For each sample create a bed of high depth regions (>150 depth to start a region and extend in both directions if >100 depth) 
- Take the UNION of all the bed files across the 10 samples (NB – only 5/10 samples were Males, so  
- For each region, annotate avgDepth and maxDepth in each of the samples 
- Filter for regions with at least 3/10 samples with MaxDepth > 250  

Additionally any homopolymer or dinucleotide repeat of >=30 bases is added to the problematic regions file. We merge repeat regions that have a gap no larger than 20 bases. If a repeat overlaps a problematic region as defined above, then we simply extend the problematic region is include the repeat region. 

Each region in the problematic regions file is further annotated with the following information 
- Overlaping genes (exonic regions only) 
- G,C,T,A count
- DistanceToTelomere
- DistanceToN 
 
In hg38, 152 genes in total have some overlap with the problematic regions file, many of which are proximate to the centromere.  3 genes of clinical cancer interest have some overlap:  FOXP1 (driver - intronic only), COL6A3 (fusion partner - intronic only) and DNMT3A (3’UTR).   All 3 cases are short PolyA regions <100 bases and will likely benefit from removing discordant pairs / soft clipped reads in these regions.  

## Known Issues / Future improvements

**Duplicate marking**
- **Read length for unpaired** - Currently this is not checked, and hence we tend to over collapse duplciates, but in some scenarios it may be valuable
- **Reads with unmapped mates** – Currently marked as duplicates based on the coordinates (and UMI) of the aligned read only.  Could lead to over-clustering.
- **Reads with mates with multiple similar local alignments** – These are currently under-clustered and lead to counting umi groups multiple times 
- **Distinguish optical vs PCR duplicates** - Duplicates should be marked as ‘optical’ if the tile distance < opticalThreshold or otherwise as PCR duplicates. 
- **Supplementary and Primary mixed up** - A fragment with a supplementary can be duplicated sometimes where there are the same 2 alignments for the read but the opposite alignment is marked as supplementary in each.  This leads us to fail to realise it is the same fragment

**UMI matching** 
- **Fixed UMI sets** - (eg TWIST) we don’t explicitly model, but doing so could lead to improvements.
- **Indel UMI errors** - will cause alignments to be different and we will fail to mark as duplicates 
- **UMI Base qual** - Not currently used, but could add value. Used in fgbio.
- **G>T errors on 1st base on UMIs** – We have observed this frequently in TWIST data, but don’t know why. 
- **INDELS near 5' read end** - We should add a final pass to UMI collapsing to allow for indels at the 5’ end of the read where we collapse groups with identical UMIs that match coordinates on 1 end and within 8 bases on the other end.
- **R2 UMI sometimes is Poly G** - Sometimes R2 fails to be read by illumin sequencers and returns Poly G.  In this case we may end up with a 2nd fragment covering R1 with an unmapped R2.

**Consensus Fragment**
- **DUAL strand consensus ref preference** - For DUAL strand consensus a qual 11 ref on one strand can cause a qual 37 alt on the other to be set to consensus ref.  This may be sub-optimal.
- **Max QUAL not always representative** - Sometimes we may see n fragments collaped with n-1 having base qual 11 and the nth fragment having base qual 37.  The low median qual indicates a potential issue, but we take the max qual right now.  This can lead to over confident DUPLEX calls. Even if we took the median instead, if one strand has just a single read with high base qual, the final consensus will have high base qual.

**Problematic regions definitions**
- REDUX should trinculeotide repeats of at least 30 length and all dinculeotide / single base repeats of  of 20-30 bases to the problematic regions file. 
- REDUX should increase minimum 10 bases outside of problematic region to 20.
- REDUX should unmap any read with discordant mate if 'repeat trimmed length' < 30 bases

 ## Version History and Download Links
