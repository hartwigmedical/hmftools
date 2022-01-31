
# Somatic Alterations in Genome (SAGE)

SAGE is a precise and highly sensitive somatic SNV, MNV and small INDEL caller.   It has been optimised for 100x tumor / 40x normal coverage, but has a flexible set of filters that can be adapted to lower depth coverage.

Key features include:
  - 4 tiered (`HOTSPOT`,`PANEL`, `HIGH_CONFIDENCE`, `LOW_CONFIDENCE`) calling allows high sensitivity calling in regions of high prior likelihood including hotspots in low mappability regions such as HIST2H3C K28M
  - kmer based model which determines a unique [read context](#read-context) for the variant + 10 bases of anchoring flanks and rigorously checks for partial or full evidence in tumor and normal regardless of local mapping alignment
  - Modified [quality score](#modified-tumor-quality-score) incorporates different sources of error (MAPQ, BASEQ, edge distance, improper pair, distance from ref genome, repeat sequencing errors) without hard cutoffs
  - Explicit modelling of ‘jitter’ sequencing errors in microsatellite allows improved sensitivity in microsatellites while ignoring common sequencing errors
  - No cutoff for homopolymer repeat length for improved INDEL handling 
  - [Phasing](#6-phasing) of somatic + somatic and somatic + germline up to ~75 bases
  - Native MNV handling 
  - Tumor sample only support
  - Multiple tumor sample support - a 'tumor' in SAGE is any sample in which we search for candidate variants and determine variant support.
  - Additional reference sample support - a 'reference' sample in SAGE is a sample in which we don't look for candidate variants, but in which we still determine variant support and read depth at each candidate location.  One potential case is to have a paired RNA sample as an additional reference to measure RNA support for candidate variants
  - Mitochondrial calling
  - An internal [alt specific base quality recalibration](#1-alt-specific-base-quality-recalibration) method
  
## Germline mode

Sage can be run in a germline mode.  See details [here](https://github.com/hartwigmedical/hmftools/blob/master/sage/GERMLINE.md).

## BAM Requirements
BAM records that are flagged as unmapped, duplicateRead or secondary/supplementary are ignored. 

Optional NM tag (edit distance to the reference) is used in the quality calculation where available otherwise it is calculated on the fly.
More information about the tag available [here](https://samtools.github.io/hts-specs/SAMtags.pdf).

While SAGE does support CRAM files, we strongly recommend converting them to BAM first as SAGE makes multiple passes over the supplied alignment files. 
Converting them first up front saves significant CPU time overall. 

## Installation

To install, download the latest compiled jar file from the [download links](#version-history-and-download-links). 

37 and 38 resources are available to download from [HMFTools-Resources > SAGE](https://resources.hartwigmedicalfoundation.nl/). 

R is used to generate the base quality recalibration charts. Required packages include `ggplot2`,`tidyr` and `dplyr`. 
R is not required if the charts are disable with the `-bqr_plot false` argument. 

# Usage

## Mandatory Arguments

Argument | Description 
---|---
tumor | Comma separated names of the tumor sample
tumor_bam | Comma separated paths to indexed tumor BAM file
out | Name of the output VCF
ref_genome | Path to reference genome fasta file
ref_genome_version | One of `37` or `38`
hotspots | Path to hotspots vcf
panel_bed | Path to panel bed
high_confidence_bed | Path to high confidence bed
ensembl_data_dir | Path to Ensembl data cache

The cardinality of `tumor` must match `tumor_bam`. At least one tumor must be supplied.

## Optional Arguments
Argument | Default | Description 
---|---|---
threads | 2 | Number of threads to use
reference | NA | Comma separated names of the reference sample
reference_bam | NA | Comma separated paths to indexed reference BAM file
chr | NA | Limit sage to comma separated list of chromosomes
max_read_depth | 1000 | Maximum number of reads to look for evidence of any `HIGH_CONFIDENCE` or `LOW_CONFIDENCE` variant. Reads in excess of this are ignored.  
max_read_depth_panel | 100,000 | Maximum number of reads to look for evidence of any `HOTSPOT` or `PANEL` variant. Reads in excess of this are ignored.  
max_realignment_depth | 1000 | Do not look for evidence of realigned variant if its read depth exceeds this value
min_map_quality | 10 | Min mapping quality to apply to non-hotspot variants
coverage_bed | NA | Write file with counts of depth of each base of the supplied bed file
validation_stringency | STRICT | SAM validation strategy: STRICT, SILENT, LENIENT

The cardinality of `reference` must match `reference_bam`.

## Optional Base Quality Recalibration Arguments

The following arguments control the [alt specific base quality recalibration](#1-alt-specific-base-quality-recalibration) logic.

Argument | Default | Description 
---|---|---
bqr_enabled | true | Enable base quality recalibration
bqr_plot | true | Plot BQR charts
bqr_sample_size | 2,000,000 | Sample size of each autosome
bqr_max_alt_count | 3 | Max support of variant before it is considered likely to be real and not a sequencing error
bqr_min_map_qual | 10 | Min mapping quality of bam record

## Optional Quality Arguments

The following arguments are used to calculate the [modified tumor quality score](#modified-tumor-quality-score)

Argument | Default | Description 
---|---|---
jitter_penalty | 0.25 | Penalty to apply to qual score when read context matches with jitter
jitter_min_repeat_count | 3 | Minimum repeat count before applying jitter penalty
base_qual_fixed_penalty | 12 | Fixed penalty to apply to base quality
map_qual_fixed_penalty | 15 | Fixed penalty to apply to map quality
map_qual_improper_pair_penalty | 15 | Penalty to apply to map qual when SAM record does not have the ProperPair flag
map_qual_read_events_penalty | 8 | Penalty to apply to map qual for additional events in read

## Example Usage

Minimum set of arguments (running in tumor only mode):

```
java -Xms4G -Xmx32G -cp sage.jar com.hartwig.hmftools.sage.SageApplication \
    -tumor COLO829v003T -tumor_bam /path/to/COLO829v003T.bam \
    -ref_genome_version 37 \
    -ref_genome /path/to/refGenome.fasta \
    -hotspots /path/to/KnownHotspots.37.vcf.gz \
    -panel_bed /path/to/ActionableCodingPanel.somatic.37.bed.gz \
    -high_confidence_bed /path/to/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed \
    -ensembl_data_dir /path_to_ensmebl_cache/ \
    -out /path/to/COLO829v003.sage.vcf.gz
```

Typical arguments running in paired tumor-normal mode:

```
java -Xms4G -Xmx32G -cp sage.jar com.hartwig.hmftools.sage.SageApplication \
    -threads 16 
    -reference COLO829v003R -reference_bam /path/to/COLO829v003R.bam \
    -tumor COLO829v003T -tumor_bam /path/to/COLO829v003T.bam \
    -ref_genome_version 37 \
    -ref_genome /path/to/refGenome.fasta \
    -hotspots /path/to/KnownHotspots.37.vcf.gz \
    -panel_bed /path/to/ActionableCodingPanel.somatic.37.bed.gz \
    -high_confidence_bed /path/to/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed \
    -ensembl_data_dir /path_to_ensmebl_cache/ \
    -out /path/to/COLO829v003.sage.vcf.gz
```

# Append Reference Samples
It is possible to append additional reference samples to an existing SAGE VCF file. A typical use case would be to append RNA without having to rerun all of SAGE.

In append mode SAGE only performs the [alt specific base quality recalibration](#1-alt-specific-base-quality-recalibration) and [normal counts and quality](#4-normal-counts-and-quality) steps.
The supplied SAGE VCF is used to determine the candidate variants and no changes are made to tumor counts, filters, phasing, de-duplication or realignment.

## Usage

## Mandatory Arguments

Argument | Description 
---|---
reference | Comma separated names of the reference sample
reference_bam | Comma separated paths to indexed reference BAM file
input_vcf | Name of the existing SAGE 2.4+ VCF
out | Name of the output VCF
ref_genome | Path to reference genome fasta file

The cardinality of `reference` must match `reference_bam` and must not already exist in the input VCF.

## Optional Arguments
Argument | Default | Description 
---|---|---
threads | 2 | Number of threads to use
chr | NA | Limit sage to comma separated list of chromosomes
max_read_depth | 1000 | Maximum number of reads to look for evidence of any `HIGH_CONFIDENCE` or `LOW_CONFIDENCE` variant. Reads in excess of this are ignored.  
max_read_depth_panel | 100,000 | Maximum number of reads to look for evidence of any `HOTSPOT` or `PANEL` variant. Reads in excess of this are ignored.  
max_realignment_depth | 1000 | Do not look for evidence of realigned variant if its read depth exceeds this value
min_map_quality | 10 | Min mapping quality to apply to non-hotspot variants

The optional [base quality recalibration](#optional-base-quality-recalibration-arguments) and [quality](#optional-quality-arguments) arguments also apply.  

## Example Usage

Minimum set of arguments:

```
java -Xms4G -Xmx32G -cp sage.jar com.hartwig.hmftools.sage.SageAppendApplication \
    -reference COLO829v003RNA -reference_bam /path/to/COLO829v003RNA.bam \
    -ref_genome /path/to/refGenome.fasta \
    -input_vcf /path/to/COLO829v003.sage.vcf.gz
    -out /path/to/COLO829v003.sage.rna.vcf.gz
```

 # Read context 
 The read context of a variant is the region surrounding it in the read where it was found.
 It must be sufficiently large to uniquely identify the variant from both the reference and other possible variants at that location regardless of local alignment.
 SAGE uses the read context to search for evidence supporting the variant and calculate the allelic depth and frequency.
 
 The core read context is a distinct set of bases surrounding a variant after accounting for any microhomology in the read and any repeats in either the read or ref genome.
 A 'repeat' in this context, is defined as having 1 - 10 bases repeated at least 2 times. 
 The core is a minimum of 5 bases long.  
 
 For a SNV/MNV in a non-repeat sequence this will just be the alternate base(s) with 2 bases either side. 
 For a SNV/MNV in a repeat, the entire repeat will be included as well as one base on either side, eg 'TAAAAAC'.
 
 A DEL always includes the bases on either side of the deleted sequence. 
 If the delete is part of a microhomology or repeat sequence, this will also be included in the core read context.
 
 An INSERT always includes the base to the left of the insert as well as the new sequence. 
 As with a DEL, the core read context will be extended to include any repeats and/or microhomology.

The complete read context is the core read context flanked on either side by an additional 10 bases. 
 
The following example illustrate how we construct and use a read context for a simple T > A SNV.  

The read context core is the variant itself expanded to cover at least 5 bases. 
Typically we use 10 bases for the flank, but for this illustration we then use an additional 5 bases on either side to get the complete read context. 
  
<pre>
Reference:                ...ACCATGGATACCATCATCACATACGA...
Variant:                                  <b>A</b>
Core read context:                      <b>CAACA</b>
Flanked read context:              <b>GATACCAACATAACA</b>
</pre>

In the following table we match the read context against bam reads in numerous ways. 
A `FULL` match includes both flanks, a `PARTIAL` match is if the read is truncated over one of the flanks but matches what is remaining, and a `CORE` match is only the core. 
A `REALIGNED` match must include both flanks but just be offset. An `ALT` match matches only the variant. All types of matches contribute to the VAF but only `FULL` and `PARTIAL` matches contribute to the `QUAL` score.

<pre>
Reference:                ...ACCATGGATACCATCATAACATACGA...
Variant:                                  <b>A</b>
Read context:                      <b>GATACCAACATAACA</b>
Full Match:               ...ACCATG<b>GATACCAACATAACA</b>TACGA...
Partial Match:                       <b>TACCAACATAACA</b>TACGA...
Core Match:               ...ACCATGGAC<b>ACCAACATAACA</b>TAACATACGA...
Realigned Match:          ...ACCCATG<b>GATACCAACATAACA</b>TACG...
Alt Match:                ...ACCATGGATACC<b>GAG</b>ATAACATACGA...
</pre>

If the variant itself is the ref (regardless of whatever else happens in the core), the read supports the ref.
 
<pre>
Reference:                ...ACCATGGATACCATCATAACATACGA...
Variant:                                  <b>A</b>
Read context:                      GATACCA<b>A</b>CATAACA
No Match:                 ...ACCATGGATACC<b>TC</b>CATAACATACGA...
No Match:                 ...ACCATGGATACCA<b>CT</b>ATAACATACGA...
Ref Match:                ...ACCATGGATACC<b>T</b>TCATAACATACGA...
Ref Match:                ...ACCATGGATACCAT<b>T</b>ATAACATACGA...
</pre>

The importance of capturing the microhomology is demonstrated in the following example. This delete of 4 bases in a AAAC microhomology is nominally left aligned as 7: AAAAC > A but can equally be represented as 8:AAACA > A, 9:AACAA > A, 10: ACAAA > A, 11: CAAAC > C etc. 
 
Using a (bolded) read context of `CAAAAACAAACAAACAAT` spanning the microhomology matches every alt but not the ref:
 
 <pre>
 REF:   GTCTCAAAAACAAACAAACAAACAATAAAAAAC 
 ALT:   GTCT<b>CAA    AAACAAACAAACAAT</b>AAAAAAC
 ALT:   GTCT<b>CAAA    AACAAACAAACAAT</b>AAAAAAC
 ALT:   GTCT<b>CAAAA    ACAAACAAACAAT</b>AAAAAAC
 ALT:   GTCT<b>CAAAAA    CAAACAAACAAT</b>AAAAAAC
 ALT:   GTCT<b>CAAAAAC    AAACAAACAAT</b>AAAAAAC
 ALT:   GTCT<b>CAAAAACA    AACAAACAAT</b>AAAAAAC
 ALT:   GTCT<b>CAAAAACAA    ACAAACAAT</b>AAAAAAC
 ALT:   GTCT<b>CAAAAACAAA    CAAACAAT</b>AAAAAAC
 ALT:   GTCT<b>CAAAAACAAAC    AAACAAT</b>AAAAAAC
 ALT:   GTCT<b>CAAAAACAAACA    AACAAT</b>AAAAAAC
 ALT:   GTCT<b>CAAAAACAAACAA    ACAAT</b>AAAAAAC
 ALT:   GTCT<b>CAAAAACAAACAAA    CAAT</b>AAAAAAC
 ALT:   GTCT<b>CAAAAACAAACAAAC    AAT</b>AAAAAAC
 ALT:   GTCT<b>CAAAAACAAACAAACA    AT</b>AAAAAAC
 ALT:   GTCT<b>CAAAAACAAACAAACAA    T</b>AAAAAAC
 </pre>
 
A similar principle applies to any repeat sequences. Spanning them in the read context permits matching with alternate alignments.
 
# Algorithm

There are 9 key steps in the SAGE algorithm described in detail below:
  1. [Alt Specific Base Quality Recalibration](#1-alt-specific-base-quality-recalibration)
  2. [Candidate Variants](#2-candidate-variants)
  3. [Tumor Counts and Quality](#3-tumor-counts-and-quality)
  4. [Normal Counts and Quality](#4-normal-counts-and-quality)
  5. [Soft Filter](#5-soft-filters)
  6. [Phasing](#6-phasing)
  7. [De-duplication](#7-de-duplication)
  8. [Realignment](#8-realignment)
  9. [Gene Panel Coverage](#9-gene-panel-coverage)

## 1. Alt Specific Base Quality Recalibration

SAGE includes a base quality recalibration method to adjust sequencer reported base qualities to empirically observed values since we observe that qualities for certain base contexts and alts can be systematically over or under estimated which can cause either false positives or poor sensitivity respectively.
This idea is inspired by the GATK BQSR tool, but instead of using a covariate model we create a direct lookup table for base quality adjustments. 
The recalibration is unique per sample.

The empirical base quality is measured in each reference and tumor sample for each {trinucleotide context, alt, sequencer reported base qual} combination and an adjustment is calculated.   This is performed by sampling a 2M base window from each autosome and counting the number of mismatches per {trinucleotide context, alt, sequencer reported base qual}.
Sites with 4 or more ALT reads are excluded from consideration as they may harbour a genuine germline or somatic variant rather than errors.    

Note that the definition of this recalibrated base quality is slightly different to the sequencer base quality, since it is the probability of making a specific ALT error given a trinucleotide sequence, whereas the sequencer base quality is the probability of making any error at the base in question.   Since the chance of making an error to a specific base is lower than the chance of making it to a random base, the ALT specific base quality will generally be higher even if the sequencer base quality matches the empirical distribution.

For all SNV and MNV calls the base quality is adjusted to the empirically observed value before determining the quality. 
SAGE produces both a file output and QC chart which show the magnitude of the base quality adjustment applied for each {trinucleotide context, alt, sequencer reported base qual} combination.
These files are written into the same directory as the output file.

A typical example of the chart is shown below. Note that each bar represents the amount that will be added to the sequencer Phred score: 

![Base Quality Adjustment](src/main/resources/readme/COLO829v003T.bqr.png)

Base quality recalibration is enabled by default but can be disabled by supplying including the`-bqr_enabled false` argument.

The base quality recalibration chart can be independently disabled by including the `-bqr_plot false` argument.
 
## 2. Candidate Variants
In this first pass of the tumor BAM(s), SAGE looks for candidate variants.
Valid candidates include a complete read context in addition to raw counts of ref and alt support (`RAD`) and their respective base quality contributions (`RABQ`).  
The raw values are calculated directly from the aligner without any filters or quality requirements.

INDELS are located using the `I` and `D` flag in the CIGAR.
SNVs and MNVs are located by comparing the bases in every aligned region (flags `M`, `X` or `=`) with the provided reference genome.
MNVs can be of up to 3 bases but with no more than one matching base between un-matching bases, ie, MNVs with CIGARs `1X1M1X` and `3X` are both considered valid MNVs of length 3.  

SAGE tallies the raw ref/alt support and base quality and selects a single read context of each variant. 
As each variant can potentially have multiple read contexts due to sequencing errors or sub-clonal populations, SAGE selects the most frequently found one as the candidate read context.
If a variant does not have at least one complete read context (including flanks) it is discarded.
All remaining variants are candidates for processing in the second pass. 

The variants at this stage have the following properties available in the VCF:

Field | Description
---|---
RC | (Core) Read Context
RC_REPS | Repeat sequence in read context
RC_REPC | Count of repeat sequence in read context
RC_MH | Microhomology in read context
RDP | Raw Depth
RAD\[0,1\] | Raw Allelic Depth \[Ref,Alt\]
RABQ\[0,1\] | Raw Allelic Base Quality \[Ref,Alt\]

Note that these raw depth values do NOT contribute to the AD, DP, QUAL or AF fields. These are calculated in the second pass. 

### Multiple Tumors
If multiple tumors are supplied, the final set of candidates is the superset of all individual tumor candidates that satisfy the hard filter criteria. 

## 3. Tumor Counts and Quality

The aim of the stage it to collect evidence of each candidate variant's read context in the tumor. 
SAGE examines every read overlapping the variant tallying matches of the read context. 
A match can be:
  - `FULL` - Core and both flanks match read at same reference location.
  - `PARTIAL` - Core and at least one flank match read fully at same position. Remaining flank matches but is truncated.  An 'N' cigar (representating a splice junction gap in RNA) may overlap both the flank and part of the core as long as the remaining flank and core match precisely.  
  - `CORE` - Core matches read but either flank doesn't.
  - `REALIGNED` - Core and both flanks match read exactly but offset from the expected position.

Failing any of the above matches, SAGE searches for matches that would occur if a repeat in the complete read context was extended or retracted.  Matches of this type we call 'jitter' and are tallied as `LENGTHENED` or `SHORTENED`. 

If the variant is not found and instead matches the ref genome at that location, the `REFERENCE` tally is incremented.

Any read which spans the core read context increments the `TOTAL` tally.  

### Modified Tumor Quality Score

If a `FULL` or `PARTIAL` match is made, we update the quality of the variant. 
No other match contributes to quality.  
There are a number of constraints to penalise the quality:
  1. as the variant approaches the edge of a read,
  2. if the read encompasses more than one variant, or
  3. if the ProperPair flag (0x02) is not set 

To do this we first calculate a modified base quality as follows:

<pre>
distanceFromReadEdge = minimum distance from either end of the complete read context to the edge of the read  
baseQuality (SNV/MNV) = BASEQ at variant location(s)  
baseQuality (Indel) = average BASEQ over core read context  
modifiedBaseQuality = min(baseQuality - `baseQualityFixedPenalty (12)` , 3 * distanceFromReadEdge) 
</pre>

We also modify the map quality taking into account the number of events, soft clipping and improper pair annotation:

<pre>
readEvents = NM tag from BAM record adjusted so that INDELs and (candidate) MNVs count as only 1 event
distanceFromReferencePenalty =  (readEvents - 1) * `map_qual_read_events_penalty (8)`^ 
improperPairPenalty = `mapQualityImproperPaidPenalty (15)`  if proper pair flag not set else 0  
modifiedMapQuality^ = MAPQ - `mapQualityFixedPenalty (15)`  - improperPairPenalty - distanceFromReferencePenalty  
</pre>

^ note that for the 6 highly polymorphic HLA genes (HLA-A,HLA-B,HLA-C,HLA-DQA1,HLA-DQB1,HLA-DQR1) we instead use modified MAPQ = min (10, MAPQ - mapQualityFixedPenalty).  We intend to improve this at some later stage by making the caller HLA type aware.


We then take the minimum of the 2 modified qualities as the read contribution to the total quality: 

<pre>
matchQuality += max(0, min(modifiedMapQuality, modifiedBaseQuality))
</pre>

A 'jitter penalty' is also calculated.  The jitter penalty is meant to model common sequencing errors whereby a repeat can be extended or contracted by 1 repeat unit.  Weakly supported variants with read contexts which differ by only 1 repeat from a true read context found in the tumor with a lot of support may be artefacts of these sequencing errors and are penalised.  If a `LENGTHENED` or `SHORTENED` jitter match is made we increment the jitter penalty as a function of the count of the repeat sequence in the microsatellite:

<pre>
`JITTER_PENALTY` += `jitterPenalty (0.25)`  * max(0, repeatCount - `jitterMinRepeatCount (3)`)
</pre>

The final quality score also takes into account jitter and is calculated as:

<pre>
`QUAL` =  matchQuality - `JITTER_PENALTY`
</pre>

### Output

The outputs of this stage are found in the VCF as:

Field | Description
---|---
RC_CNT\[0,1,2,3,4,5,6\] | Read Context Count \[`FULL`, `PARTIAL`, `CORE`, `REALIGNED`, `ALT`, `REFERENCE`, `TOTAL`\]
RC_QUAL\[0,1,2,3,4,5,6\] | Read Context Quality \[`FULL`, `PARTIAL`, `CORE`, `REALIGNED`, `ALT`,`REFERENCE`, `TOTAL`\]
RC_JIT\[0,1,2\] | Read Context Jitter \[`SHORTENED`, `LENGTHENED`, `JITTER_PENALTY`\]
AD\[0,1\] | Allelic Depth  (=\[RC_CNT\[5\], RC_CNT\[0\] + RC_CNT\[1\] + RC_CNT\[2\] + RC_CNT\[3\] + RC_CNT\[4\]\] )
DP | Read Depth (=RC_CNT\[6\])
AF | Allelic Frequency (=AD\[1\] / DP)
QUAL | Variant Quality (=RC_QUAL\[0\] + RC_QUAL\[1\] - RC_JIT\[2\])

### Hard Filters

To reduce processing the following hard filters are applied: 

Filter | Default Value | Field
---|---|---
hard_min_tumor_qual |30| `QUAL`
hard_min_tumor_raw_alt_support |2| `RAD[1]`
hard_min_tumor_raw_base_quality |0| `RABQ[1]`
filtered_max_normal_alt_support |3| Normal `AD[1]`

Note that hotspots are never hard-filtered.

The first 3 filters are excluded from this point onwards and have no further processing applied to them.  
The filtered_max_normal_alt_support is applied at the final step of the algorithm, solely to reduce file size. 
The filtered_max_normal_alt_support does not apply to germline variants in the same local phase set as passing somatic variants.
 
## 4. Normal Counts and Quality

Evidence of each candidate variant is collected in all of the supplied reference bams in the same manner as step 2. 

RNA bams are valid reference sources.

## 5. Soft Filters

Given evidence of the variants in the tumor and normal we apply somatic filters. 
The key principles behind the filters are ensuring sufficient support for the variant (minimum VAF and score) in the tumor sample and validating that the variant is highly unlikely to be present in the normal sample.

The filters are tiered to maximise sensitivity in regions of high prior likelihood for variants. 
A hotspot panel of 10,000 specific variants are set to the highest sensitivity (TIER=`HOTSPOT`) followed by medium sensitivity for exonic and splice regions for the canonical transcripts of a panel of cancer related genes (TIER =`PANEL`) and more aggressive filtering genome wide in both high confidence (TIER=`HIGH_CONFIDENCE`) and low confidence (TIER=`LOW_CONFIDENCE`) regions to ensure a low false positive rate genome wide.   These tiers can be customised by providing alternative bed files as configuration

The specific filters and default settings for each tier are:

Filter  | Hotspot | Panel | High Confidence | Low Confidence | Field
---|---|---|---|---|---
min_tumor_qual**|70***|100|160|240|`QUAL`
min_tumor_vaf|0.5%|1.5%|2.5%|2.5%|`AF`
min_germline_depth|0|0|10 | 10 | Normal `RC_CNT[6]`
min_germline_depth_allosome|0|0|6 | 6 | Normal `RC_CNT[6]`
max_germline_vaf****|10%|4%|4% | 4% | Normal`RC_CNT[0+1+2+3+4]` / `RC_CNT[6]`
max_germline_rel_raw_base_qual|50%|4%|4% | 4% | Normal `RABQ[1]` / Tumor `RABQ[1]` 

** These min_tumor_qual cutoffs should be set lower for lower depth samples.  For example for 30x tumor coverage, we recommend (Hotspot=40;Panel=60;HC=100;LC=150)

*** Even if tumor qual score cutoff is not met, hotspots are also called so long as tumor vaf >= 0.08 and  allelic depth in tumor supporting the ALT >= 8 reads.  This allows calling of pathogenic hotspots even in known poor mappability regions, eg. HIST2H3C K28M.

**** A special filter (max_germline_alt_support) is applied for MNVs such that it is filtered if 1 or more read in the germline contains evidence of the variant.

If multiple tumors are supplied, a variant remains unfiltered if it is unfiltered for any single tumor. 

The germline criteria are only evaluated against the primary reference, ie, the first in the supplied reference list.
If no reference bams supplied, the germline criteria are not evaluated.

Soft filters can be disabled using the `disable_soft_filter` parameter.

Soft filters become hard filters when the `hard_filter` flag is included. 

To set the parameters at the command line append the tier to the filter eg `hotspot_min_tumor_qual` and `high_confidence_min_tumor_qual` set the value of the min_tumor_qual for the `HOTSPOT` and `HIGH_CONFIDENCE` tiers respectively.

## 6. Phasing

### Local Phase Set
Local phasing implies that two or more variants co-exist on the same read. 
In practice, somatic variants are phased using the read contexts of nearby germline or other somatic variants. 
As we are comparing read contexts (with an additional buffer of 50 bases) we are only able to phase variants within a maximum of 60 bases of each other.
This distance is further reduced if the read contexts we are comparing have been extended in opposing directions due to repeats or microhomology.

Phasing is interesting for somatic calling from 3 perspectives: 
  - understanding the somatic mutational mechanism which has led to the variant; and 
  - understanding the functional impact of the variation.
  - prediction of neoeptiopes for MHC class 1 presentation
  
Regarding mechanism, multiple somatic cis-phased variants can frequently occur together with prominent mechanisms being 2 base MNVs (eg. CC>TT in UV Signatures and CC>AA in Lung smoking signatures) and micro-inversions (which are normally called as offsetting INS and DEL).

Phasing somatic and germline variants together can also aid in understanding somatic mechanisms - for example if the germline has a 6 base deletion in a microsatellite and the tumor has a 7 base deletion, then the likely somatic mechanism is a 1 base deletion. 

Phasing is also important for functional impact particularly in 2 prominent cases: 
  - 2 nearby phased frameshift variants can lead to an inframe (and potentially activating) INDEL; and 
  - 2 phased synonymous SNVs in the same codon could potentially cause a nonsense or missense effect together.    

Two variants are considered phased if their read context cores and any intervening bases are identical after adjusting for their relative position.
This is demonstrated in the example below where two SNVs share an identical sequence of bases.

<pre>
REF: CAACAATCGAACGATATAAATCTGAAA
A>T: CAACAATCGA<b>T</b>CGATACAATC
T>C:       TCGATCGATA<b>C</b>AAATCTGAAA
</pre>

Similarly, SNVs, MNVs and INDELs may be phased together. Any variants that are phased together are given a shared `LPS` (local phase set) identifier.

If multiple tumors are supplied, phasing is evaluated only on the primary tumor, ie, the first in the supplied tumor list.

### Local Realignment Set
Local realignment implies that two (or more) variants are equivalent. 
As the aligner assesses each read independently, small changes in the position of a variant within a read can lead to different interpretations.
The following example illustrates this:

<pre>
POS: 123456789...
REF: AAATGATTT...
ALT: AAAT<b>A</b>ATTT...
</pre>

This is shown above as SNV 5:G>A but can equally be represented as the following phased INDEL 4:TG>T and SNV 7:T>A so long as the subsequent bases match:

<pre>
POS: 123456789...
REF: AAATGATTT...
ALT: AAAT A<b>A</b>T...
</pre>

We can detect local realigned variants using a similar process to phasing but without adjusting the relative position by the INDEL insert/delete sequence.

Any variants that are can be locally realigned are given a shared `LRS` (local realigned set) identifier.

### Phased Inframe Indels

If two phased frameshift variant in a single coding exon together form an inframe INDEL, then both are given a shared `PII` (phased inframe indel) identifier.

## 7. De-duplication

### Realigned Indels

While the read context is designed to capture a unique sequence of bases, it it sometimes possible that repeat sequences in the flanks of the read context coupled with an aligners alternate view on the same event can cause duplicate calls of the same event.    If SAGE finds two phased INDELs of the same type at the same position where one is a subset of the other, then the longer is filtered with `dedup`.

Alternative alignments near read edges can also cause INDELs to be mistaken for SNV.  For example, in the local realignment set example above, the read could be represented by a either a single SNV, or a phased INDEL and SNV.  In such cases, we keep the variant with the highest individual quality as well as any variants also phased with it and filter any others with `dedup`.

### SNV / MNV

Any passing SNVs that are phased with and contribute to a passing somatic MNV are filtered with `dedup`. 

If the MNV is comprised of both somatic and germline SNVs, additional logic applies. 
Typically, the MNV will be filtered, except when the MNV is in a coding region and impacts more than one base of the same codon impacted by the SNV. 
In this case the SNV is filtered to capture the functional impact of the mixed somatic and germline variant. 

The following example illustrated this as the combined impact of the germline and somatic SNVs is different to just the somatic SNV. 
In this case we keep the MNV and filter the somatic SNV. The germline SNV remains filtered. 

Chromosome | Position | Ref | Alt | Type | Protein Impact
---|---|---|---|---|---
19 | 43,382,367 | GT | AG | Mixed | p.Thr43Leu
19 | 43,382,367 | G | A | Somatic | p.Thr43Ile
19 | 43,382,368 | T | G | Germline | p.Thr43Pro

If the MNV is comprised of only germline SNVs but does not appear itself at all in the germline, it remains unfiltered.  
 
Any MNVs that have a germline component and all associated SNVs (including somatic) are given a shared `MSG` (mixed somatic germline) identifier.



## 8. Realignment

Inframe indels with microhomology are re-aligned to the right if the left-aligned variant is not in a coding region but the right-aligned variant is.

For example, because of the `AG` microhomology at this (v37) location, the following KIT variants are equivalent but the first will be interpreted as a splice variant while the second will be interpreted as an inframe missense.

```
4:55593579 CAGAAACCCATGTATGAAGTACAGTGGA > C
4:55593581 GAAACCCATGTATGAAGTACAGTGGAAG > G
```

Similarly, the following EGFR inserts are equivalent:

```
7:55248980 C > CTCCAGGAAGCCT
7:55248992 T > TTCCAGGAAGCCT
```


## 9. Gene Panel Coverage
To provide confidence that there is sufficient depth in the gene panel a count of depth of each base in the gene panel is calculated and written to file for each tumor sample. 
This can be disabled by setting the `panel_coverage` parameter to `false`.

The file shows the number of bases with 0 to 30 reads and then buckets reads in intervals of 10 up to 100+.

For a sample with approximately 30x depth this may appear as: 

```
gene	0	1	2	...	27	28	29	30-39	40-49	50-59	60-69	70-79	80-89	90-99	100+
BRCA1	0	0	0	...	23	54	116	1854	2834	875	5	0	0	0	0
ERBB2	0	0	0	...	106	154	203	2315	832	83	0	0	0	0	0
TP53	0	0	0	...	31	41	59	636	192	35	0	0	0	0	0
```

While a 100x depth sample might appear something more like:

```
gene	0	1	2	...	27	28	29	30-39	40-49	50-59	60-69	70-79	80-89	90-99	100+
BRCA1	0	0	0	...	0	0	0	0	0	0	0	8	390	1700	3672
ERBB2	0	0	0	...	0	0	0	0	0	18	257	590	1311	1111	611
TP53	0	0	0	...	0	0	0	0	0	0	0	90	423	343	376
```

A 'missed variant likelihood' is calculated using poisson as the mean probability of not finding at least 3 reads coverage for an allele given the depth distribution over the whole gene.

# Variant Pipeline
A number of post processing steps are applied to the SAGE output.

## PON Filtering
To eliminate recurrent variants and artifacts we constructed a Panel of Normal (PON) by first running SAGE over 1000 germline samples and recording any variants with at least 3 reads and sum of base quality > 30. 
The frequency (`PON_COUNT`) of each variant and the maximum read support in any one sample ('PON_MAX') are then aggregated into the PON file.  
Variants with more than 1 observation are retained.

We use the PON file to filter SAGE output of any variant that appears in more than 6 samples for the LOW_CONFIDENCE & HIGH_CONFIDENCE & PANEL tiers and 10 samples for HOTSPOTS.  
Additionally, PANEL and HOTSPOT variants must have at least 1 sample with 5 or more reads support to be PON filtered.

The annotation and filtering are done with BCFTools using commands similar to the following:

```
bcftools annotate -a SageGermlinePon.1000x.37.vcf.gz -c PON_COUNT,PON_MAX $sage_vcf -O z -o $annotated_vcf

bcftools filter -e 'PON_COUNT!="." && INFO/TIER="HOTSPOT" && PON_MAX>=5 && PON_COUNT >= 10' -s PON -m+ $annotated_vcf -O u | \
bcftools filter -e 'PON_COUNT!="." && INFO/TIER="PANEL" && PON_MAX>=5 && PON_COUNT >= 6' -s PON -m+ -O u | \
bcftools filter -e 'PON_COUNT!="." && INFO/TIER!="HOTSPOT" && INFO/TIER!="PANEL" && PON_COUNT >= 6' -s PON -m+ -O z -o $pon_filtered_vcf
```

As the 38 PON is smaller with only 98 germline samples it is recommended that the following filters are used:

```
bcftools annotate -a SageGermlinePon.98x.38.vcf.gz -c PON_COUNT,PON_MAX $sage_vcf -O z -o $annotated_vcf

bcftools filter -e 'PON_COUNT!="." && INFO/TIER="HOTSPOT" && PON_MAX>=5 && PON_COUNT >= 5' -s PON -m+ $annotated_vcf -O u | \
bcftools filter -e 'PON_COUNT!="." && INFO/TIER="PANEL" && PON_MAX>=5 && PON_COUNT >= 2' -s PON -m+ -O u | \
bcftools filter -e 'PON_COUNT!="." && INFO/TIER!="HOTSPOT" && INFO/TIER!="PANEL" && PON_COUNT >= 2' -s PON -m+ -O z -o $pon_filtered_vcf
```

## SnpEff
SnpEff is run over the output and then summarised for each variant with a separate post processing application.

# Performance Characteristics
Performance numbers were taken from a 72 core machine using paired normal tumor COLO829 data with an average read depth of 35 and 93 in the normal and tumor respectively. 
Elapsed time is measured in minutes. 
CPU time is minutes spent in user mode. 
Peak memory is measured in gigabytes.

Threads | Elapsed Time| CPU Time | Peak Mem
---|---|---|---
1 | 696 | 751 | 10
8 | 98 | 776 | 13
16 | 62 | 873 | 13 
24 | 49 | 880 | 14 
32 | 45 | 943 | 15

# Version History and Download Links
- Upcoming
  - Renamed `assembly` to `ref_genome_version` with valid values [37, 38]
- [2.8](https://github.com/hartwigmedical/hmftools/releases/tag/sage-v2.8)
  - Right align inserts that would otherwise be outside a coding region in the same manner as deletes
- [2.7](https://github.com/hartwigmedical/hmftools/releases/tag/sage-v2.7)
  - Calculate NM field if not present in alignment file
- [2.6](https://github.com/hartwigmedical/hmftools/releases/tag/sage-v2.6)
  - Coverage is now calculated on supplied bed file rather than on panel bed file
  - Added validation_stringency parameter
- [2.5](https://github.com/hartwigmedical/hmftools/releases/tag/sage-v2.5)
  - Gene panel coverage including estimate of missed variant likelihood 
  - Changed default value of max_germline_rel_raw_base_qual for HOTSPOTS to 50% (from 100%)
  - Improved VAF estimation by including alt match type 
- [2.4](https://github.com/hartwigmedical/hmftools/releases/tag/sage-v2.4)
  - Added SageAppendApplication to [append additional reference samples](#append-reference-samples) to existing SAGE output. 
  - Do not hard filter germline variants in the same local phase set as passing somatic variants
  - Large skipped reference sections (representating a splice junction gap in RNA) contribute towards FULL or PARTIAL matches.
- [2.3](https://github.com/hartwigmedical/hmftools/releases/tag/sage-v2.3)
  - Extend local phase set detection to maximum of 60 bases
  - Favour reads with variants closer to the centre when determining read context
  - Fix bug creating BQR plots with too many quality scores
- [2.2](https://github.com/hartwigmedical/hmftools/releases/tag/sage-v2.2)
  - Realignment of inframe indels
  - Improved MNV deduplication
  - Detection of phased inframe indels
  - Base Quality Recalibration
  - Improved sensitivity in high depth regions
  - Tumor only support
  - Mitochondria support
  - Multiple tumor support
  - Multiple reference (or RNA) support
  - Removed explicit RNA support (can use additional reference instead)
  - Performance and memory improvements
- [2.1](https://github.com/hartwigmedical/hmftools/releases/tag/sage-v2.1)
  - Reduced memory footprint
  - Add version info to VCF
  - RNA support
  - CRAM support
  - Filter variants with ref containing bases other then G,A,T,C
  - Look for read context of variants that are partially soft clipped
  - Ref genome 38 support
- 2.0
  - Revamped small indel / SNV caller
