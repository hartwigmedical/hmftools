# Esvee - Structural Variant Calling

## Overview and algorithm


## BAM Prep Pre-Assembly Filtering

Prep generates a maximally filtered SV BAM file by identifying candidate SV junctions and extracting all reads that may provide support to 
that junction.

### Command

```
java -cp esvee.jar com.hartwig.hmftools.esvee.prep.PrepApplication 
  -sample 'REF_SAMPLE_ID,TUMOR_SAMPLE_ID'
  -bam_file '/sample_data/REF_SAMPLE_ID.bam,/sample_data/TUMOR_SAMPLE_ID.bam'
  -ref_genome /path_to_ref_genome_fasta/
  -ref_genome_version 38 
  -known_fusion_bed /ref_data/known_fusions.38.bedpe
  -bamtool /tools/sambamba/sambamba 
  -output_dir /sample_data/output/ 
  -threads 16
```

#### Mandatory Arguments

Argument | Description 
---|---
sample | Sample IDs separated by ','
bam_file | BAM file paths separated by ','
ref_genome | Reference genome fasta file
ref_genome_version | 37 (default) or 38
bamtool | Sambamba or Samtools, required to sort and index the output BAMs
output_dir | Output directory
threads | Thread count

#### Optional Arguments

Argument | Description 
---|---
known_fusion_bed | BED file with known fusion pair coordinates, require only 1 fragment for junctions
blacklist_bed | See below for explanation


### Algorithm

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

The blacklist is the combination of the existing encode blacklisted regions and all regions with 200x or greater depth found in 4 out of 11 reference samples (40x mean coverage) used to identify artefacts (and 4 out of 8 for HG38).
Regions that overlap PANEL or fusion KB genes are excluded from the blacklist unless the max coverage is >2000x. These regions mainly capture long repeat sections of the genome with poorly aligned reads and make up 13M bases of the genome (0.4%).


## Assembly Building



### Command

```
java -cp esvee.jar com.hartwig.hmftools.esvee.assembly.AssemblyApplication 
  -tumor TUMOR_SAMPLE_ID 
  -reference REF_SAMPLE_ID
  -tumor_bam /sample_data/TUMOR_SAMPLE_ID.bam
  -reference_bam /sample_data/REF_SAMPLE_ID.bam
  -junction_file /sample_data/output/TUMOR_SAMPLE_ID.esvee.prep.junction.tsv
  -ref_genome /path_to_ref_genome_fasta/
  -ref_genome_version 38
  -write_types 'JUNC_ASSEMBLY;PHASED_ASSEMBLY;ALIGNMENTS;BREAKEND;VCF'
  -output_dir /sample_data/output/ 
  -threads 16
```

#### Mandatory Arguments

Argument | Description 
---|---
tumor | Tumor sample ID
tumor_bam | Path to Prep tumor BAM file
reference | Tumor sample ID (can be omitted in tumo-only mode)
reference_bam | Path to Prep reference BAM file
junction_file | Path to Prep junction TSV file, assumes named as 'TUMOR_SAMPLE_ID.esvee.prep.junction.tsv'
ref_genome | Reference genome fasta file
ref_genome_version | 37 (default) or 38
write_types | Minimum required is VCF for latter steps
output_dir | Output directory
threads | Thread count

#### Optional Arguments

Argument | Description 
---|---
decoy_genome | Decoy fastq sequences file, eg use HG38 decoys for a GRCH37 run


### Algorithm



## Reference Depth Annotation

Once unfiltered variants have been identified, they are annotated with the depth matching the reference genome in each input BAM. 
This feeds into the VAF calculations in the caller routine below.

```
java -cp esvee.jar com.hartwig.hmftools.esvee.depth.DepthAnnotator \
  -sample 'REF_SAMPLE_ID,TUMOR_SAMPLE_ID'
  -bam_file '/sample_data/REF_SAMPLE_ID.bam,/sample_data/TUMOR_SAMPLE_ID.bam'
  -input_vcf TUMOR_SAMPLE_ID.esee.raw.vcf.gz
  -output_vcf TUMOR_SAMPLE_ID.esee.ref_depth.vcf.gz
  -ref_genome /path_to_ref_genome_fasta/
  -ref_genome_version 38
  -threads 16
```

#### Arguments

Argument | Description 
---|---
sample | Sample IDs separated by ','
bam_file | BAM file paths separated by ','
ref_genome | Reference genome fasta file
ref_genome_version | 37 (default) or 38
input_vcf | Input VCF from assembly, assumes named as 'TUMOR_SAMPLE_ID.esvee.raw.vcf.gz'
output_vcf | Output VCF, default is TUMOR_SAMPLE_ID.esvee.ref_depth.vcf.gz
threads | Thread count



## Variant Calling and Filtering

The final step is to filter and annotate all variants, and then to write out 3 VCFs:
- somatic VCF - SAMPLE_ID.esvee.somatic.vcf.gz
- germline VCF - SAMPLE_ID.esvee.germline.vcf.gz
- unfiltered VCF - SAMPLE_ID.esvee.unfiltered.vcf.gz


### Command

```
java -cp esvee.jar com.hartwig.hmftools.esvee.caller.CallerApplication 
  -sample TUMOR_SAMPLE_ID
  -reference REF_SAMPLE_ID
  -input_vcf /sample_data/output/TUMOR_SAMPLE_ID.esvee.ref_depth.vcf.gz
  -ref_genome_version 38
  -pon_sgl_file /ref_data/sgl_pon.38.bed.gz
  -pon_sv_file /ref_data/sv_pon.38.bedpe.gz
  -known_hotspot_file /ref_data/known_fusions.38.bedpe
  -repeat_mask_file /ref_data/repeat_mask_data.37.fa.gz
  -output_dir /sample_data/output/ 
```

#### Mandatory Arguments

Argument | Description 
---|---
sample | Tumor sample ID
reference | Reference sample ID
ref_genome_version | 37 (default) or 38
output_dir | Output directory

#### Optional Arguments

Argument | Description 
---|---
input_vcf | VCF from reference depth routine above, assumed named as 'TUMOR_SAMPLE_ID.esvee.ref_depth.vcf.gz'
pon_sgl_file | PON for SGL breakends
pon_sv_file | PON for SVs
known_hotspot_file | Known hotspot SVs, matches known-pair fusions as used by Linx
repeat_mask_file | Repeat mask file



## Known issues and future improvements


### Version History and Download Links
- [1.0](https://github.com/hartwigmedical/hmftools/releases/tag/esvee-v1.0)
