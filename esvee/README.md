# Esvee - Structural Variant Calling

## Overview and algorithm


## ESVEE Prep - pre filtering

SV Prep generates a maximally filtered SV BAM file by identifying candidate SV junctions and extracting all reads that may provide support to 
that junction. The BAM file is intended to be fed into the GRIDSS assembly.   SV Prep reduces the overall runtime of GRIDSS SV calling by ~80%.

In tumor-normal mode, SV Prep may be run first on the tumor and then a 2nd time on the reference sample using the junctions found in the tumor mode, 
to ensure all potential evidence in the reference sample is collected for candidate tumor junctions. 

Example usage of SV Prep can be found [here](https://github.com/hartwigmedical/hmftools/blob/master/pipeline/scripts/run_sv_calling)

### Running SvPrep

#### Mandatory Arguments

Argument | Description 
---|---
sample | Sample ID
bam_file | Input BAM file
ref_genome | Reference genome fasta file
ref_genome_version | 37 (default) or 38
output_dir | Output directory
threads | Thread count

#### Optional Arguments

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
java -cp esvee.jar com.hartwig.hmftools.esvee.prep.PrepApplication 
  -sample SAMPLE_ID
  -bam_file /sample_data/SAMPLE_ID.bam
  -ref_genome /path_to_ref_genome_fasta/
  -ref_genome_version [37 or 38] 
  -output_dir /sample_data/output/ 
```

The bed file resources can be downloaded from [HMFTools-Resources > DNA Pipeline](https://console.cloud.google.com/storage/browser/hmf-public/HMFtools-Resources/dna_pipeline/):
- know fusion pair locations BED (see resources 'sv')
- blacklist locations BED (see resources 'sv')

### Overview and algorithm

#### Extract Fragments

Ignoring reads that are duplicate (or where the primary is a duplicate), secondary or contain soft clipping with more than 16 consecutive 
bases of PolyG/C), and excluding sites in blacklisted regions, parse through BAM end to end to identify breakend sites with CREDIBLE soft clipping

- At least 1 read with 50 base alignment and  abs(Insert size - M length) > 5 bases (ie test for short fragments with adapter) AND a soft clip length of 30 with >=75% of soft clip bases with qual > 25. 
- Soft clip length of <30 is allowed where the soft clip meet the polyA LINE criteria (ie 16 of first 18 soft clip bases must be A/T). 
- The read must also not be a repeat expansion (ie the first 9 bases of soft clip and the last 9 bases of aligned read are not matching 1,2 or 3 nucleotide repeats).       
- At least 1 additional read which have any length soft clipping at the same base OR within 50 bases with <=1 high quality mismatch between the soft clip locations  (not required for HOTSPOT regions)   
- At least 1 read supporting directly or indirectly with a MAPQ of 20 (not required for HOTSPOT regions)

Sites with aligned INDELs of >= 32 bases are also treated as candidate sites.

Additionally, to ensure we capture SV with long (inexact) homology we also identify sites where there are >=5 reads within a 500 base region 
with insert size outside the max(1000,99.75%) insert size (or >=3 reads if insert size is more than twice that length) with their mates also 
starting within a 1000 base region of each other. The range can be further extended to the inner side of the fragment up to the 99.75% fragment 
size if more reads can be found in that region with long insert sizes and mates within 1kb of each other.  If such a candidate region is found 
and the reads do not support a site of credible soft clipping, then create a junction regardless on the innermost base of the reads supporting 
the potential breakend.

For each site we also obtain the following reads (excluding reads with alignments that overlap blacklist regions or which have PolyG/C tails) and their mates
- All reads with soft clipping that matches the orientation and position of the variant (+/-50 bases)
- All reads within min[99.75%,1kb] range fragment length of the site on the correct side of the breakend with read facing the breakend and mate is unmapped, 
interchromosomal, has the same orientation or has an insert size outside the percentile range [0.25,99.75]
- All reads that overlap the breakend and contain an INDEL of 5+ bases

For the last 2 categories (discordant & indel containing reads), we filter if they have more than max(5,25% of soft clip length) soft clip bases with base qual < 25, since these frequently cause FP calls in GRIDSS. 

A BAM is written for each input sample which contains the reads described above and their mates. This is then fed into the Esvee assembly routine described below.

## Reference Depth Annotation

Esvee Prep also has an additional feature to replace the depth annotation of GRIDSS (ie the annotation of REF and REFPAIR) with a faster implementation.  This can be run with the following command: 

```
java -cp esvee.jar com.hartwig.hmftools.esvee.depth.DepthAnnotator \
  -input_vcf ${gridss_vcf} \
  -output_vcf ${final_vcf} \
  -samples "${reference_id},${tumor_id}" \
  -bam_files "${reference_bam},${tumor_bam}" \
  -ref_genome ${ref_genome} \
  -ref_genome_version ${ref_genome_version} \
  -threads ${threads} \
```

Please see the example [script](https://github.com/hartwigmedical/hmftools/blob/master/pipeline/wgs_scripts/run_gridss) for how to run this post GRIDSS.

### SV Prep Blacklist

The blacklist is the combination of the existing encode blacklisted regions and all regions with 200x or greater depth found in 4 out of 11 reference samples (40x mean coverage) used to identify artefacts (and 4 out of 8 for HG38).   Regions that overlap PANEL or fusion KB genes are excluded from the blacklist unless the max coverage is >2000x. These  regions mainly capture long repeat sections of the genome with poorly aligned reads and make up 13M bases of the genome (0.4%).  

## Known issues and future improvements

- **Soft Clip end bias** - There are several artefacts that occurs solely on the 3' end of both R1 and R2.   We could potentially explicitly filter sites with a lot of 3' support but no 5' support.
- **Short INDEL Artefacts** - Short INDELS (10-40 bases) may be called as SGL and not filtered (normally with poor qual or shortish assemblies).  This could be addressed by identifying these better in GRIPSS
- **Microsatellites** - Related, other artefacts may still be called immediately adjacent to microsatellites.   Additional filtering may help
- **MT chromosome** - currently dropped
- **Variants near blacklisted regions** - GRIDSS will ignore any read that overlaps a blacklisted region.   Hence, we cannot call any breakpoint which is within ~30-70 bases of a blacklisted region.
- **Max softclip overlap** - Currently must be 30 bases of soft clip overlap.  This may cause us occasionally to miss 1 or 2 reads support at one end of a break junction if it is only identified at one end.

# Version History and Download Links
- [1.1](https://github.com/hartwigmedical/hmftools/releases/tag/sv-prep-v1.1)
- [1.0](https://github.com/hartwigmedical/hmftools/releases/tag/sv-prep-v1.0.1)
