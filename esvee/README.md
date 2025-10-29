# Esvee - Structural Variant Calling

## Overview

Esvee is a structural variant caller optimised for short read sequencing that identifies somatic and germline somatic rearrangements.

Esvee is run in 4 steps
- BAM Filtering & Prep
- Assembly & Alignment
- Reference Depth Annotation
- Variant Calling & Filtering 

The full algorithm for each step is described in the algorithm section below.

### Command

Esvee can be run as a single command or by running the 4 steps described below in turn. 
A single command example is:

```
java -jar esvee.jar  
  -tumor TUMOR_SAMPLE_ID 
  -reference REF_SAMPLE_ID
  -tumor_bam /sample_data/TUMOR_SAMPLE_ID.bam
  -reference_bam /sample_data/REF_SAMPLE_ID.bam
  -ref_genome /path_to_ref_genome_fasta/
  -ref_genome_version 38
  -write_types 'PREP_STANDARD;ASSEMBLY_STANDARD' 
  -known_hotspot_file /ref_data/known_fusions.38.bedpe
  -pon_sgl_file /ref_data/sgl_pon.38.bed.gz
  -pon_sv_file /ref_data/sv_pon.38.bedpe.gz
  -target_regions_bed /ref_data/CoverageCodingPanel.38.bed
  -artefact_pon_sgl_file /ref_data/sgl_artefact_pon.38.bed.gz
  -artefact_pon_sv_file /ref_data/sgl_artefact_pon.38.bed.gz
  -repeat_mask_file /ref_data/repeat_mask_data.38.fa.gz
  -bamtool /tools/sambamba/sambamba 
  -output_dir /sample_data/output/ 
  -threads 16
```


## STEP 1: BAM Filtering & Prep

Prep generates a maximally filtered SV BAM file by identifying candidate SV junctions and extracting all reads that may provide support to 
those junctions.

### Command

```
java -cp esvee.jar com.hartwig.hmftools.esvee.prep.PrepApplication 
  -sample 'TUMOR_SAMPLE_ID,REF_SAMPLE_ID'
  -bam_file '/sample_data/TUMOR_SAMPLE_ID.bam,/sample_data/REF_SAMPLE_ID.bam'
  -ref_genome /path_to_ref_genome_fasta/
  -ref_genome_version 38 
  -known_hotspot_file /ref_data/known_fusions.38.bedpe
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
bam_stringency | SAM validation strategy: STRICT, SILENT, LENIENT

## STEP 2: Assembly & Alignment

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
reference | Tumor sample ID (can be omitted in tumor-only mode)
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

## STEP 3: Reference Depth Annotation

Once unfiltered variants have been identified, they are annotated with the depth matching the reference genome in each input BAM. 
This feeds into the VAF calculations in the caller routine below.

```
java -cp esvee.jar com.hartwig.hmftools.esvee.depth.DepthAnnotator \
  -sample 'TUMOR_SAMPLE_ID,REF_SAMPLE_ID'
  -bam_file '/sample_data/TUMOR_SAMPLE_ID.bam,/sample_data/REF_SAMPLE_ID.bam'
  -input_vcf TUMOR_SAMPLE_ID.esvee.raw.vcf.gz
  -output_vcf TUMOR_SAMPLE_ID.esvee.ref_depth.vcf.gz
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

#### Optional Arguments

Argument | Description 
---|---
bam_stringency | SAM validation strategy: STRICT, SILENT, LENIENT
unmap_regions | BED file of unmapped regions in Redux

## STEP4: Variant Calling and Filtering

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
  -target_regions_bed /ref_data/CoverageCodingPanel.38.bed
  -artefact_pon_sgl_file /ref_data/sgl_artefact_pon.38.bed.gz
  -artefact_pon_sv_file /ref_data/sgl_artefact_pon.38.bed.gz
  -known_hotspot_file /ref_data/known_fusions.38.bedpe
  -repeat_mask_file /ref_data/repeat_mask_data.38.fa.gz
  -output_dir /sample_data/output/ 
```

#### Mandatory Arguments

Argument | Description 
---|---
sample | Tumor sample ID
reference | Reference sample ID
ref_genome_version | 37 (default) or 38
output_dir | Output directory
prep_dir | Directory with Prep output files

#### Optional Arguments

Argument | Description 
---|---
input_vcf | VCF from reference depth routine above, assumed named as 'TUMOR_SAMPLE_ID.esvee.ref_depth.vcf.gz'
pon_sgl_file | PON for SGL breakends
pon_sv_file | PON for SVs
target_regions_bed | Path to BED file for targeted panel
artefact_pon_sgl_file | Artefact PON for SGL breakends
artefact_pon_sv_file | Artefact PON for SVs
known_hotspot_file | Known hotspot SVs, matches known-pair fusions as used by Linx
repeat_mask_file | Repeat mask file

## Algorithm

### Key concepts 
Esvee is a structural variant caller that evaluates evidence for candidate breakends in a BAM and outputs forms fully assembled, aligned junctions, fully annotated and filtered based on read support.  
 
Some key terms and concepts in Esvee are set out below: 

Term | Definition 
---|---
Breakend | A candidate position and orientation for a structural variant 
Split read | A read which directly overlaps a breakend and is normally soft clipped at or near the breakend (depending on homology) 
Discordant read | A read with mate which is unmapped OR on different chromosome OR Fragment Length > 99.75% or < 0.25% of distribution. Discordant reads <1kb from and facing a junction are kept by SV_PREP and considered in assembly extension. 
INDEL read | A read with an INDEL that exceeds the minLength threshold of SV (32 bases). These are used to form and support junctions.  Shorter INDELs can support a split breakend (i.e. one formed from split reads) if offering an equivalent per-base representation
Junction | A linked pair of breakends that together form a structural variant 
Chained breakends | A pair of breakends that face each other and are explicitly cis-phased by fragments. This normally represents a short-templated piece of DNA which may be inserted into another location. 
Phase group | Maximal set of breakends that share at least 2 fragments with another breakend in the group. Breakends which may imply short DEL/DUP/INS may also be added to the same phase group, even without shared fragment evidence. Phase groups are created to limit complexity of assembly: all assembly merging and extension is done within phase groups only. 
Remote region | Co-ordinates and orientation of a set of overlapping (or near for gaps) reads that are paired with assembled or discordant reads. Breakends may have one or many remote regions. 
LINE insertion site | A site where a mobile element has been inserted characterised by 2 opposite orientation breakends either facing each other with less than 30 bases gap or overlapping with less than 40 bases and with either a long PolyA insertion sequence on the negatively oriented breakend or a long polyT insertion sequence on the positively oriented breakend. The insertion sequence size may vary from just a few bases to a long sequence and is often unmappable due to the repetitive mobile elements in the human genome.  These insertion sites have special logic throughout each of the Esvee steps to maximise sensitivity. 
SV Length | For INV, DEL, DUP and INS, the SV length is defined as the position difference between the 2 locations (excluding the last mapped base) + length of any insert sequence.  

### Key thresholds and tolerances
Much of the logic depends in Esvee depends on assembling reads into contiguous sequence, first locally and then merging assemblies and remote regions within phase groups. The algorithms use the following tolerances throughout Esvee:
- **Minimum length**: We require 32 bases to call a variant. At the Esvee-Prep and local breakend assembly stage we also require a soft clip of at least 32 bases to retain a junction (an exception is made for regions with high discordant fragment support or LINE insertion sites) 
- **Low Quality SNV errors**: Low quality mismatches (rawBQ<26) are ignored and deemed to always match an existing assembly 
- **Minimum assembly overlap**: We require 50 bases of overlap to merge and extend assemblies (with an exact matching 20 base seed) OR 20 bases to merge and extend reference bases
- **Mismatch tolerance**: When comparing reads and assemblies we allow 0 high quality mismatches for sequences of < 15 bases, 1 high quality mismatches for sequences of 15 to 100 bases, and then 1 additional high mismatches for each additional 200 bases of sequence overlap more than 100 bases. Note that a 1 or 2 base indel counts as 1 mismatch. A difference in unit count in a microsatellite repeat (1 count diff for repeats <= 10, 2 count diff for repeats 11-15, 3 count diff for repeats 16+) also count as 1 mismatch. When extending reference bases we always allow up to 2 mismatches

Esvee also employees concepts of a modified alignment scores and MAPQs to try to represent absolute versus relative likelihood of mismatch.  THese are defined as follows:

#### Adjust alignment score (AdjAS) 
Modified length of an alignment after allowing for inexact homology, repeats and mismatches. Defined as: 
```
adjAS = AS - inexact homology length – repeatBases[repeatCount>2]
```
Note that the alignment score incorporates both length and number of mismatches and that the repeats are only evaluated for the bases outside of the homology region (so that this is not double counted) 
#### Modified MAPQ (modMAPQ)  
The modified MAPQ is intended to convert the MAPQ which is a relative MAPQ into something more akin to an ‘absolute’ MAPQ.  The MAPQ is penalised if the length is short or the alignment score is low relative to the length.
```
modMAPQ = MAPQ * min [1, adjAS/max(100,alignLength)]^2 
```

### STEP 1: BAM Filtering & Prep 
This step generates a set of maximally filtered SV BAM files and an initial set of candidate SV junctions from input tumor/normal BAM files. The SV BAM files includes all candidate split and discordant reads that are proximate to candidate junctions that may provide support to that junction. 

Ignoring reads with 50%+ low qual bases (<26), are duplicate (or where the primary is a duplicate), secondary or contain soft clipping with more than 16 consecutive bases of PolyG/C), Prep parses through the BAMS end to end to identify breakend sites with CREDIBLE soft clipping: 
- At least 1 read with AdjAlignmentScore > 40 and abs(Insert size - M length) > 5 bases (ie test for short fragments with adapter) AND a soft clip length of >32 with >=75% of soft clip bases with qual > 25. We also tolerate AS < 40 if AS >= 35 and AS – XS >= 15.
- Soft clip length of <32 is allowed where the soft clip meet the polyA LINE criteria (ie 16 of first 18 soft clip bases must be A/T). 
- The read must also not be a simple repeat expansion. ie the first 9 bases of soft clip and the last 9 bases of aligned read are not matching 1,2 or 3 nucleotide repeats. 
- At least 1 additional read which have soft clipping of any length at the same base OR within 50 bases with <=1 high quality mismatch between the soft clip locations (not required for HOTSPOT regions)
  
Sites with aligned INDELs of >= 32 bases are also treated as candidate sites. 

Additionally, to ensure we capture SV with long (inexact) homology we also identify sites where there are >=5 reads within a 500 base region with insert size outside the max(1000,99.75%) insert size (or >=3 reads if insert size is more than twice that length) with their mates also starting within a 1000 base region of each other. The sites must also be on the same chromosome within +/1MB with a DEL or DUP orientation or create a link between a known hotspot region.  At least one read must have MAPQ>=40 and adjAS>75. The range can be further extended to the inner side of the fragment up to the 99.75% fragment size if more reads can be found in that region with long insert sizes and mates within 1kb of each other. If such a candidate region is found and the reads do not support a site of credible soft clipping, then create a breakend regardless on the innermost base of the reads supporting the potential breakend. 

Sites with an AF below 0.5% (1% for discordant) are filtered.

For each site we also obtain the following reads and their mates:
- All reads with soft clipping that matches the orientation and position of the variant (+/-50 bases) 
- All reads within 99.75% range fragment length of the site on the correct side of the breakend with read facing the breakend and mate is unmapped, inter-chromosomal, has the same orientation or has an insert size outside the percentile range [0.25,99.75] 
- All reads that overlap the breakend and contain an INDEL of 3+ bases 
For the last 2 categories (discordant & indel containing reads), we filter if they have more than max(5,25% of soft clip length) soft clip bases with base qual <=25, since these frequently cause FP calls in GRIDSS. 

A BAM is written for each input sample which contains the reads described above and their mates.  

Esvee also counts discordant fragments genome wide in various buckets and writes to a separate counts file. These categories are:
* Translocations
* Inversions (<1kb, 1-5kb, 5-100kb, 100+kb)
* Deletions (1-5kb, 5-100kb, 100+kb)
* Duplications (1-5kb, 5-100kb, 100+kb)

### STEP 2A: Local breakend assembly  
The candidate breakends generated by Esvee prep are organised into breakend groups, covering a set of breakends within 1kb of each other. All reads from the Esvee-Prep BAMs are retrieved for each breakend group.  

Each read has the following adjustments made to it prior to assembly: 
- **Drop low qual reads** - any reads with >50% of low qual bases (<26) are dropped altogether 
- **PolyG trimming** – any continuous stretch of 4 or more Gs (Cs on reverse-strand reads) are trimmed from the read. These bases and quals are entirely truncated from the read. 
- **Quality trimming** – Trim low quality 3' softclip tails such that any 3' base tail taken from the trimmed read has at least 65% of bases at qual 26+. To maximise LINE insertion sensitivity, the first 18 bases of the SC are NOT trimmed regardless of QUAL IF at least 16 of them are either A or T AND the bases are not an extension of an aligned PolyA/T sequence of at least 8 bases in the reference. 
- **Indel to soft-clips** – any INDEL of length >= 3 and < minBases [32] allows a read to set indel-inferred unclipped start and end read positions, to allow these reads to offer support for split junctions.

Discordant fragments with unmapped mates are ignored in subsequent steps if unmapped trimmed read length < 50 or if MAPQ = 0.  
#### Assembly 
Breakend supporting reads identified by Prep are used to build initial breakend assemblies. These must either be soft clipped on the junction side and cross the breakend or indel containing reads with indels >= min length bases (currently 32). Discordant reads are not included at this stage as we have not yet confirmed their mates are part of the same assembly. 

In the first step the ‘extension’ sequence is built from the breakend outwards. For a sequence built from softclipping, we use the consensus of soft-clipped bases from reads with clipping at the correct orientation and position. For a sequence built from an indel, we use the consensus of reads containing the relevant indel cigar element at the same position.

These initial assemblies are filtered unless at least 1 read has a soft clip > 32 bases and at least one other read has a soft clip >= 16 bases OR the sequence contains a PolyA/PolyT LINE sequence of 16 bases and a 2nd read with supporting soft clip > = 8 bases with an insertion site orientation. These two key anchoring reads also must have different 5' fragment coordinates. However, note that candidate junctions with LINE source insertion orientations are filtered as they may be promiscuous causing issues downstream. These are always assembled from the LINE insertion site end.

Each of the remaining candidate overlapping supporting reads is now tested against these unique breakend-extension assemblies. These reads will contribute to initial junction support counts, but not be used to form the initial extension sequence. For softclip junctions this includes indel containing reads, reads with a softclip of the correct orientation within 2bp, and reads with SNVs past the junction location rather than clipping. For indel junctions this includes reads with clipping past the junction, or any indel causing the read's indel-inferred unclipped start/end to cross the junction.

At this stage, if enough reads differ from the consensus with high qual mismatches exceeding the mismatch tolerance the extension can be potentially split into two junctions (specifically, if each junction has at least max(20% of total, 5) supporting reads).

To allow calling of variants with long homology which may not have soft clipping, if no assembly can be created with sufficiently long soft clip, but there exists a remote region with more than 3 discordant fragments supporting the link then an assembly is made using only the fragments associated with that specific remote region. Only reads that support the specific remote location are allowed to form the assembly  The breakend is arbitrarily initially placed 32 bases from the innermost base of the assembly. Only reads with MAPQ>20 are considered in the assembly construction. 

All remaining breakend assemblies are then extended into the reference base direction from their supporting reads. Only junction split reads with 10+ soft-clipped bases and at least 2 (4 if junction flanks an adjacent ref microsatellite) high quality soft clipped bases matching the consensus, can be used to build the initial assembly ref base sequence. 

#### Assembly deduplication 
Duplicate assemblies may exist at either the same breakend or within +/-50 bases with the same orientation. Assemblies are deduped (the assembly with the higher initial read support is retained) if they satisfy an assembly comparison within mismatch tolerances (see 'Key thresholds and tolerances'). We will prioritise keeping split over discordant junctions, and prioritise keeping junctions with germline support where possible.

#### Local assembly filtering 
Any local breakend which fully aligns to a decoy or ALT (excluding HLA ALTs) sequence (top alignment with <=3 mismatches and AS/length > 0.9) in the hg38 genome are filtered at this stage and not processed further.  

For any junction with average of support reads less than MAPQ < 10, the reference base sequence is aligned.  If the aligner returns MAPQ < 10 and no XA tag we filter the assembly as it cannot lead to a breakend with our later alignment rules. 

### STEP 2B: Phasing 
Phase groups are created by maximally linking any breakends which share at least 2 fragments. Since short DEL, DUP and INS will not share discordant reads on either side of the junction, if there are no supplementary reads that directly support the junction, they will not share any reads.  Hence, we also merge any 2 proximate breakends into the same phase group if: 
- they form a DUP orientation <500b OR DEL orientation < 1kb AND 
- both breakends have at least 1 split read with concordant mate on the soft clipped side OR at least one side has a PolyA / PolyT tail sequence with insertion site orientation.

During phasing, all candidate remote linking sites are collected for each breakend. These are taken from discordant reads, breakend assembly read mates and breakend read supplementaries, and are established from the remote read or supplementaries coordinates (ie its chromosome and read start and end alignments). We only allow remote region creation from candidate discordant fragments if the local non junction read has AS >= 75 and lowQualPerc < 25%. Remote reads with overlapping alignments are merged into sets of remote regions, and then cached against each assembly and are used later in phasing and assembly merging.

At this stage each phase group consists of multiple breakends:
- Sharing fragment support, or  
- Are proximate to each other (for short DEL, DUP, INS) 
- These breakends may be single breakends, breakend pairs, or represent chained/complex structural variants. 

Each breakend consists of: 
- A local assembly and the reads that make up that assembly 
- A set of local discordant reads proximate to that breakend which may support the breakend 
- A set of remote regions consisting of the mates of the assembled and local discordant reads 
- A set of unmapped mates of the assembled reads and local discordant reads. 
- A phase group which contains any breakends which are proximate or with shared fragment support 
 
### STEP 2C: Merge & extend assemblies 
Breakend assemblies are merged and extended to form junctions. There are 6 sub-steps: 
<img width="536" alt="image" src="https://github.com/user-attachments/assets/a680046e-5ba0-4e79-9c47-6c1bf592a6ff">

In each step, Esvee attempts to merge existing assemblies with other assemblies, remote regions and/or mates of assembled reads. Assemblies are merged if a suitable 20 base exact seed can be matched, and the seed can be extended to the end of each assembly sequence with the minimum required overlap and up to the prescribed number of high- quality mismatches (see 'Key thresholds and tolerances'). When 2 breakends are assembled, if there is a gap between the reference sequences, an insert sequence is recorded. If there is an overlap of reference sequences, the reference bases are also extended.

#### 1. Assemble local pairs 
Since short DEL & DUP are the most common variant types, Esvee first tries to merge proximate breakends (within 500 bases for DUP or 1000 bases for DEL) with consistent orientations that may be a DEL or DUP, prioritising pairs that share indel containing reads.

#### 2. Local alignments 
Each unlinked breakend is next checked if it can be resolved by local alignment. The same merge criteria is applied, but instead of aligning between 2 local assemblies, Esvee instead aligns the breakend assembly to the local reference genome within +/- 500 bases. If an alignment is made and the variant length > minLength a new breakend is created at the alignment location. If an alignment is made and the variant length < minLength, the candidate breakend is filtered.

#### 3. Assemble distal junctions 
Within a phase group, Esvee next attempts to merge any pair of unlinked distal breakend assemblies OR individual breakends with discordant and unmapped mates, that may form a paired junction assembly and are directly linked by same fragment support. Before attempting this match, the extension sequences are extended with reads from remote and unmapped regions as follows:
* Check for reads that have a 20bp exact match to the existing sequence, and a 50bp match with usual mismatch tolerances from 'Key thresholds and tolerances'. We allow down to a 30bp total match length if the sequence is non-repetitive
* If one or more matching reads exist, perform an internal assembly consensus between them and eliminate reads failing usual mismatch tolerances
* Use the remaining reads to extend the extension sequence
* If extension was performed, attempt this process again (as reads that didn't match the existing extension sequence may match now)

After doing this, the breakend pairs are merged based on matching the sequence assembly. When this occurs, the assembly is also maximally extended on both ends from the junction outwards, including in the reference sequence with local junction mates if the the required overlap is found. After all distal mergers are made, remaining candidate breakends are tested to see if they can be merged into and extend the existing merged assemblies.

Breakends which are merged into already formed junctions are marked as ‘SECONDARY’. Secondary links are not attempted for junctions already established as local indels unless they have a LINE tail. Also, assemblies are branched into multiple assemblies if there is at least 5 reads and 20% maximum support supporting an alternative alignment. This may occur in the case of foldback inversions.

For LINE insertions sites with 2 breakends, both breakends are assembled simultaneously.

#### 4. Chained assemblies 
Esvee next merges chained or facing breakends that share fragments and are at least 30 bases apart and less than 1kb. To merge, soft clipping in the reference base must match the facing breakend, and the supporting fragments must be concordant along the phased assembly.
 
#### 5. Infer local gaps 
If any unassembled mates remain for any unlinked assemblies which are within 200 bases of an assembly on the 3’ side of the read and do not have softclipping in the direction of the assembly junction, then extend the assembly to include the read sequence and assume the bases in between match the reference sequence.

#### 6. Collapse duplicates 
Finally we do two rounds of deduplication:
* We drop any previously duplicated branched assemblies that never formed a chain
* We merge two phased assemblies into a single assembly if their sequences can be matched within tolerances. This is only allowed if either: 
- The matching sequence overlaps at least 1 junction on either assembly OR 
- The 2 sequences are either end of the same LINE insertion site

Ahead of alignment, we also filter any still-unlinked assembly that is likely to just be a local indel below our min variant length. Specifically, these are:
* Softclip-based assembly, not linked to anything else in its phase group
* Extension length < 85bp
* More than 1/3 of assembly support is from indel reads, OR the two supporting reads with longest assembly length are both indel reades

### STEP 2D: Alignment & variant calling 

#### Alignment 
Esvee now has a set of unique assemblies which may relate to a single candidate breakend, a junction pair or a complex set of chained breakends. Each unique assembly is aligned using BWA-mem, and the gap open and mismatch penalties are also raised from 4 to 6 to incentivise producing softclips where edges of an alignment are very marginal.

BWA may return one primary alignment as well as one or more supplementary alignments. Since BWA can assign an unreliable MAPQ to supplementary alignments, any supplementary alignments are realigned again using BWA with the primary alignment of the re-query kept and any further supplementaries dropped. Similarly, since BWA can be disinclined to create a supplementary alignment for a long softclip (50bp+) that is part of a very long single alignment, we manually requery the softclip portion in these cases.

At this stage a modified MAPQ is calculated for each alignment 
```
modMAPQ = MAPQ * min [1, AdjustedAlignmentScore/max(100,alignLength - IHOM length)]^2
```
where:
```
AdjustedAlignmentScore = Alignment score – IHOM length – repeatBases[repeatCount>2] 
```
 
Note that if (alignmentScore + 15 < 0.77 * (length – inexact homology length)) the modMAPQ is set to 0. This helps to filter long but biologically implausible alignments. If the alignment contains a single indel at or above the min variant length (32), then the alignment score penalty associated with this indel is recovered for the purpose of this test, as well as in the adjusted alignment score. Two other additional criteria apply:
* if the assembly is for a cigar-based indel and BWA has returned two alignments, then we tolerate the alignments without checking modMAPQ as long as raw MAPQ >= 5
* otherwise, if an alignment has AdjustedAlignmentScore < 50, then reject it unless has NM < 4 and is within 1kb of a 50+ length alignment

The interpretation of the alignment depends on both the modified map quality and the XA tag which will display the alternative alignments if there are a small number of alternatives. Assemblies with no alignments or with all alignments with modMAPQ < 10 and NULL XA tags are ignored.

For assemblies with a single alignment, a junction is inferred from CIGAR and added to the VCF if the alignment has modMAPQ>=10 and the alignment has an ‘I’ or ‘D’ CIGAR element of 32+ bases with at least 50 bases of matched bases (M cigar) on either side after trimming. A single breakend is created in the VCF if there is either a soft clip length of at least 32 bases or a PolyA/PolyT INS sequence of at least 1.5x length of PolyA in the adjacent reference bases), and either the modMAPQ>=10 OR the alignment has a non NULL XA tag.

For assemblies with 2+ alignments, for each pair of consecutive alignments in the assembly with modMAPQ >=10, Esvee writes a pair of breakends representing a junction to a VCF. For alignments with modMAPQ < 10 AND XA!=NULL, the default and alternative alignments are checked to see if they can make a short variant (length < 100kb) with adjacent alignments, and if so the alignment which makes the shortest variant is chosen. If there are consecutive alignments with ModMAPQ<10 and XA!=NULL, adjacent alignments with modMAPQ >=10 are checked first, followed by any remaining consecutive alignments with XA!=NULL. If 2 identical length variants are possible between consecutive alignments with modMAPQ < 10 and XA!=NULL, then prefer the default alignment. If no short variant can be made but either the XA!=NULL or MAPQ>=5 then use the default alignment. Otherwise, the alignment is converted to insert sequence. A summary of the behaviour is shown in the below table: 
 
Case | Alignment behaviour  
---|---
Has XATag with <100kb alignment | Use shortest (or default if equal). Add 15 to QUAL 
Has XATag and MAPQ > 0 No short alignment | Use default alignment 
No XA Tag. modMAPQ>=5 | Use default alignment 
Other | No alignment (convert to insert sequence) 

For each alignment with modMAPQ<10 and either (XATag!=NULL OR alignment is in DUX4 region) regardless of whether the alignment was inferred or was converted to an insert sequence, then all the alternatives are recorded in the INSALN field for any adjacent breakends or the ALTALN field if an alignment.  
 
Note that 2 special exceptions are made for clinically relevant fusions in low mappability: 
- **SSX2**: SS18 Exon 10 to SSX2 Exon 6 mutations are common pathogenic fusions in Synovial Sarcomas, but may be confused with the homolog SSX2B.  Therefore, for any alignment which falls in the range of intron 5 of SSX2 (hg19: X:52,729,628-52,731,680; hg38 chrX: 52,700,578-52,702,630) and SSX2B (hg19: X: 52,784,877-52,786,929; hg38 chrX:52,755,800-52,757,852) with a downstream genic orientation with MAPQ < 20 and the assembly has a breakend mapped to another chromosome is assumed to have the SSX2 alignment. 
- **DUX4**: DUX4 is a special case as it has known clinically relevant fusions but many identical copies in the reference genome. DUX4 regions are defined as {GL000228.1:20000-125000, 4:190930000-191030000, 10:135420000-135520000} on hg19 and {chr4:190060000-190190000, chr10:133660000-133770000} in hg38. If there are 2 or more consecutive modMAPQ<3 alignments, then the INSALN is merged for those cases.

A final filtering pass is done to eliminate spurious alignments resulting from likely single read extensions. If all the following conditions are satisfied we drop all breakends corresponding to the phased assembly:
- Phased assembly consists of only a initial assembly, without links
- Assembly is softclip based and does not contain a poly-A/T LINE tail
- The alignment of this phased assembly does not resolve as a short del or dup (<=1kb)
- The longest extension read is >=50bp, but the second longest is <32bp
 
#### Homology annotation and precise breakend alignment 

Where alignments are overlapping, the junction is aligned to the midpoint of the subset of overlapping bases that require that lead to the least total number of mismatches to the reference genome across both breakends. If there are multiple non-contiguous regions with identical number of mismatches, then the longest is selected. The confidence interval CIPOS is set to the bounds of the lowest mismatch region and the bases in this region are set to the homology sequence (HOMSEQ). The bounds of the overlap relative to the aligned breakpoint is set to be the inexact homology confidence interval (IHOMPOS). In the case of DUPs called off alignments with I cigar elements, the leftmost valid junction position is chosen instead, and inexact homology is set based on exact homology.

#### QUAL Annotation 

The QUAL is set to the SUM of the alignment qualities at either end, with an adjustment to penalise alignments with low read count support: 
```
QUAL = [modMAPQ(localBE) + modMAPQ(remoteBE)]*[readSupport / (readSupport + halfQualSupport)] 
```
HalfQualSupport is set to 4 by default. 

If a breakend is in multiple assemblies, the qual is calculated across each of the assemblies and the highest QUAL is used. 

#### Read count annotations 
For each breakend the total variant fragments (VF) that overlap the breakend across all assemblies is recorded. This is also broken down into fragments with reads that directly overlap the breakend (SR) and discordant fragments that span but do not overlap the breakend. Note that spanning fragments are NOT included in counts for INS or DEL and DUP < 1000 bases in length as we cannot uniquely determine if they are discordant. 

### STEP 3: Reference counts & filtering 
#### Reference count annotations 
For each breakend, Esvee queries the full BAM files and annotates the fragments with a read that directly overlaps the breakpoint (REFSR) and spans the breakpoint (REFRP). These are also used to calculate a variant allele frequency (AF) for the variant. To be consistent, for INS or DEL and DUP < 1000 bases in length as well as SGLs, the REFRP is also ignored from the denominator of the AF calculation. If the VAF is > 90% and the variant falls into an ‘unmapped region’ as defined in REDUX, then the REF_DEPTH is assumed to be the sample median. 

### STEP 4: Filtering 
The following filters are applied to the variant with a context of ‘any sample’ (ie PASS if any sample meets the criteria) or ‘all samples’ (pass if the criteria is met across all samples for the variant: 

Filter Name | Samples | Definition | Junction | LINE Site | Single | Hotspot
---|---|---|---|---|---|---
minQual | Any | QUAL  | 30 (WGS) / 60 (PANEL) |  30 (WGS) / 60 (PANEL) <sup>1</sup> |  30 (WGS) / 60 (PANEL)  | 30 
minSupport  | Any   | VF | 4 | 4<sup>1</sup>  | 6 | 2 
minAF | Any | min(AF[BE1],AF[BE2]) | 0.001 | 0.001 | 0.05  | 0.001 
invShortLowVafHom<sup>5 | Any | min(AF[BE1],AF[BE2])  | 3<=IHOMLEN<6:  min(0.1,200 * shortINVRate); IHOMLEN>=6 min(0.2,400 * shortINVRate) <sup>5</sup> | NA | NA  | NA 
invShortFragLowVaf<sup>7 | Any | min(AF[BE1],AF[BE2])  | min(0.05, 50 * shortFragmentINVRate) | NA | NA | NA
delShortLowVaf<sup>8 | Any | min(AF[BE1],AF[BE2])  | 0.05 | NA | NA | NA 
minLength<sup>2</sup>  | All | EndPos-StartPos+InsSeqLength | 32 | NA | NA | 32
shortFrags | All | Lengthmedian - NumSD * LengthstdDev/sqrt(VF)<sup>3</sup>   | 3 | NA | NA | 3 
minAnchorLength | All | AlignLength – repeatLength – Homology | 50 | NA | 50<sup>4</sup>  | 50 
sbArtefact<sup>6</sup> | All | SB | NA | NA | 1.0 | NA
invShortIsolated | All (targeted panel only) | Candidate inversion distance to nearest non-artefact breakend<sup>9</sup> | 100 | NA | NA | NA

<sup>1. For pairs of SGL breakends which resemble a likely LINE insertion site (see above) the SUM(Qual) is used for both breakends. </sup>

<sup>2. Same chromosome junctions only. </sup>

<sup>3. implies the sampled average fragment length should be within 3 standard deviations of the sample median length (note the cutoff is also capped at 0.6*SD below median length).  Standard deviation is estimated as Lengthmedian-length16th percentile </sup>

<sup>4. For SGL BE, the inserted sequence length must also meet these requirements </sup>

<sup>5. Only applied to variants with type=INV and LEN<3kb. ShortINVRate = proportion of fragments genome wide that support a short INV < 5kb </sup>

<sup>6. Only for SGL 5' end contains GCCGTATCATTAAAAA or GTAGATCTCGGTGGTC OR SGL 3' containing TTTTTAATGATACGGC or GACCACCGAGATCTAC </sup>

<sup>7. Only applied to variants with type=INV and LEN<300bp. ShortINVRate = proportion of fragments genome wide that support a short INV < 1kb </sup>

<sup>8. Only applied to variants with type=DEL, IHOMLEN>=10, AVGLEN+LEN >= 1.5*99.75th percentile of fragment lengths, and LEN<3kb </sup>

<sup>9. Candidate inversions are inversions <100bp that are unchained or have unique fragment positions (UFP) < 4. Eligible non-artefact breakends face and are within 100bp of the candidate inversion's breakend, and are not filtered by PON or any inversion-specific filters </sup>

Note that for pairs of breakends at LINE insertion sites, if one has a PASS filter we should always PASS the other side.   

In targeted mode, SGL breakends are only retained if in a targeted region. Paired-end variants are only retained if either breakend is within 1kb of a targeted region.

#### Variant Deduplication 
If the same precise breakend is found to PASS multiple times in the VCF then retain the variant with the highest QUAL only 

#### Germline or Somatic determination 
A consolidated VCF is produced showing all soft filters. If a germline sample is present and BOTH max(germline AF/TumorAF) > 0.1 AND germlineAD / tumorAD > 0.01 the variant is deemed to be germline, else somatic. Separate vcfs are written for PASS and PON somatic and germline variants only (in tumor only mode just a somatic vcf filter is written). A PON filter is also applied to the somatic variant vcf based on both a SGL BE PON (for SGL breakends) and a paired BE PON (for junctions). SGL breakends with candidate or alterate alignments that fall in the SSX2, SSX2B or DUX4 regions are never PON filtered. For pairs of breakends at LINE insertion sites, if one variant is marked as germline, then both should be considered as germline.  

## Summary of LINE insertion site behaviour

Stage  | Special rules 
---|---
Esvee Prep| min 32 base length is not required 
Trimming | Don’t trim the first 18 bases if at least 16 of them are PolyA or T  
Local assembly | Prioritise reads with the 5’ in the softclip and which reach beyond the PolyA for extension; Allow any length of PolyA to match. Set the length in the assembly to be the median of lengths with additional bases beyond the PolyA or else just the longest PolyA sequence if none exist; Require only 16 bases longest and 8 bases 2nd longest to retain soft clip 
Phasing | Phase breakends even if neither has a locally concordant mate. Also use a permissive matching routine to merge the assemblies
Assembly extension | Allow remote regions and unmapped mates of both sides to extend the assembly 
Alignment | Call single breakend if PolyA length exceeds 1.5x reference PolyA length (no 50 base minimum)
Filters | MinLength, ShortFrags, minAnchorLength filters not applied; PASS if either side PASSes ; MinSupport & MinQual uses qual of both sides 
Germline vs Somatic | Mark as germline if either side meets germline filters 

## Output 

### VCF INFO fields
Field |Description
---|---
ALTALN |Potential alternative alignments of segment in the format chr:start|strand|cigar|mapq
ASMID	|Unique id(s) of assembly(s) containing the breakend
ASMLEN	|Total length(s) of assembly(s) containing the breakend
ASMLNKS |Breakend id of breakends linked by assembly
ASMSEG	| #(s) of segments in assembly(s) containing the breakend
AVGLEN	|Average implied length of fragments supporting the junction
BEAOR	|Breakend orientation(s) in assembly(s)
BEAPOS	|Breakend position(s) in assembly(s)
BEOR	|Breakend orientation(s) in reference genome
CIPOS	|Confidence interval around breakend position (for homology)
HOMSEQ	|Homology sequence at junction
HOTSPOT	|Is Known fusion hotspot
IHOMPOS	|Offset positions of inexact homology
INSALN	|Potential alignment locations of insert sequence in the format chr:start|strand|cigar|mapq. Populated when max MAPQ =0 and <= 5 alternative alignments
INSRMP	|Portion of inserted sequence whose alignment overlaps the repeatmasker repeat
INSRMRC	|IInserted sequence repeatmasker repeat class
INSRMRT	|Inserted sequence repeatmasker repeat type
LINE	|LINE Insertion Site
MATEID	|Id of other breakend in junction
PON_COUNT	|PON count if in PON
SEGALEN	|Aligned length of segment(s) in reference genome
SEGID	|Unique id(s) of segment(s) containing the breakend
SEGMAPQ	|MAPQ of segment containing the breakend with highest QUAL contribution
SEGRL	|Repeat length of segment  with highest QUAL contribution
SEGSCO	|Alignment score of segments containing the breakend with highest QUAL contribution
SVID	|ID shared by both breakends in the variant
SVTYPE	|Type of structural variant
UFP |Unique fragment positions
MLR| Max local repeat

### VCF sample specific fields
Field |Description
---|---
AD	| Allelic depths for the ref and alt alleles in the order listed
AF	|	Allele frequency of the breakend
DF	|	Count of discordant fragments with a read either side of the breakend
DP	|	Approximate read depth
GT	|	Genotype
REF	|	Count of fragments supporting the reference with a read overlapping the breakend
REFPAIR	|	Count of paired fragments supporting the ref with a read either side of the breakend
SB	|	Proportion of split reads with 3' end facing the breakend. Fragments with both reads split are counted in both directions
SF	|	Count of fragments supporting the breakend with a read overlapping the breakend
VF	|	Total variant fragments supporting the breakend

## Known issues and future improvements
Know sources of errors
- Misintepretation of INDELs in long repeats can sometimes cause poor quality consensus sequences in assembly extension.
- SGL AF will be systematically underestimated if we cannot extend the assembly.
- **Somatically activated LINE insertion sites** - Some insertion sites of LINE elements may become active LINE source elements themselves. These may appear to be BOTH insertion and source sites for LINE elements and may lead to overcounting of support at insertion sites. 

Alignment
- We should analyse additional supplementary alignments arising from re-query of initial supplementary alignments (currently dropping)

Esvee has some implicit and explicit assumptions on reads, qualities and alignments:
- **AS field** - AS is currently required
- **Low qual masking** -  assumes a high proportion of bases have qual > 30 
- **Read lengths** - Soft clip & alignment score assumptions require read lengths > 80 bases
- **Fragment lengths** - We use 1000,500 to refer to short DEL, DUP respectively.  Ideally this should depend on fragment lengths.
- **Low Qual INDELS** - Technologies with many low quality indel errors (eg Ultima) may have assembly impacted. These should be masked from assembly
- **Hard clipping** - Esvee prep may not retain reads with hard clipping at or near junctions

Other planned improvements
- Variant visualisations
- Multi tumor and reference sample support including donor support
- Append mode
- Investigate region specific filtering for IG,TCR & HLA regions
- Switch minSupport from MAX to SUM of samples & lower minSupport to 3. Could analyse given recent duplicate collapsing improvements in Redux and Esvee prep

### GRIDSS and GRIPSS SV filtering
- [Old Gripss information](./README_GRIPSS.md)

### Version History and Download Links
- [1.0](https://github.com/hartwigmedical/hmftools/releases/tag/esvee-v1.0)
