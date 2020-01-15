
# Somatic Alterations in Genome (SAGE)
SAGE is a somatic SNV, MNV and small INDEL caller with features including:
  - kmer based model which determines a unique [read context](#read-context) for the variant + 25 bases of anchoring flanks and rigorously checks for partial or full evidence in tumor and normal regardless of local mapping alignment
  - modified [quality score](#quality) incorporates different sources of error (MAPQ, BASEQ, edge distance, improper pair, distance from ref genome, repeat sequencing errors) without hard cutoffs
  - Explicit modelling of ‘jitter’ sequencing errors in microsatellite allows improved sensitivity in microsatelites while ignoring common sequencing errors
  - no cutoff for homopolymer repeat length for improved INDEL handling 
  - 3 tiered (Hotspot,Panel,Wide) calling allows high sensitivity calling in regions of high prior likelihood including hotspots in low mappability regions such as HIST2H3C K28M
  - [Phasing](#6.-phasing) of somatic + somatic and somatic + germline up to 25 bases
  - Native [MNV handling](#mnv-handling) 

 # Read context 
 
 The core read context is the distinct set of bases surrounding a variant after accounting for any microhomology in the read and any repeats in the read or ref genome.
 In this context, a repeat is defined as having 1 - 10 bases repeated at least 2 times. 
 The core is a minimum of 5 bases long.  
 
 For a SNV in a non-repeat sequence this will just be the single alternate base with 2 bases either side. 
 For a SNV in a repeat, the entire repeat will be included as well as one base on either side, eg 'TAAAAC'.
 
 A DEL always includes the bases on either side of the deleted sequence. 
 If the delete is part of a microhomology or repeat sequence, this will also be included in the core read context.
 
 An INSERT always includes the base to the left of the insert as well as the new sequence. 
 As with a DEL, the core read context will be extended to include any repeats and/or microhomology.
 
 The importance of capturing the microhomology is demonstrated in the following example. This delete of 4 bases in a AAAC microhomology is  
 nominally left aligned as 7: AAAAC > A but can equally be represented as 8:AAACA > A, 9:AACAA > A, 10: ACAAA > A, 11: CAAAC > C etc. 
 
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
 
 A similar principle applies to any repeat sequences. Spanning them in the read context permits matching alternate alignments.
 
 The complete read context is the core read context flanked on either side by an additional 25 bases. 
 
# Algorithm

There are 8 key steps in the SAGE algorithm described in detail below:
  1. [Candidate Variants And Read Contexts](#1.-candidate-variants-and-read-contexts)
  2. [Tumor Counts and Quality](#2.-tumor-counts-and-quality)
  3. [Hard Filters](#3.-hard-filters)
  4. [Normal Counts and Quality](#4.-normal-counts-and-quality)
  5. [Soft Filter](#5.-soft-filters)
  6. [Phasing](#6.-phasing)
  7. [MNV Handling](#7.-mnv-handling)
  8. [Output](#8.-output)
 
 
## 1. Candidate Variants And Read Contexts

In this first parse of the tumor BAM, SAGE uses the `I` and `D` flag in the CIGAR to find INDELs and compares the bases in every aligned 
region (flags `M`, `X` or `=`) with the provided reference genome to find SNVs. 
Note that there are no base quality or mapping quality requirements when looking for candidates.

SAGE tallies the raw ref/alt support and base quality and collects the read contexts of each variant.
Once finished, each variant is assigned its most frequently found read context as its primary one. 
If a variant does not have at least one complete read context (including flanks) it is discarded.
All remaining variants are then considered candidates for processing in the second pass. 

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

Note that these values do NOT contribute to the AD, DP, QUAL or AF fields. These are calculated in the second pass. 

## 2. Tumor Counts and Quality

The aim of the stage it to collect evidence of each candidate variant's read context in the tumor. 
SAGE examines every read overlapping the variant tallying matches of the read context. 
A match can be:
  - `FULL` - Core and both flanks match read at same reference location.
  - `PARTIAL` - Core and at least one flank match read fully at same position. Remaining flank matches but is truncated. 
  - `CORE` - Core matches read but neither flank does.
  - `REALIGNED` - Core and both flanks match read exactly but offset from the expected position.

Failing any of the above matches, SAGE searches for matches that would occur if a microsatellite in the complete read context was extended or retracted. 
Matches of this type we call jitter and are tallied as `LENGTHENED` or `SHORTENED`.  

Lastly, if the base the variant location matches the ref genome, the `REFERENCE` tally is incremented while any read which spans the core read context increments the `TOTAL` tally. 

### Quality

If a `FULL` or `PARTIAL` match is made, we update the quality of the variant. 
No other match contributes to quality.  
There are a number of constraints to penalise the quality if it:
  1. approaches the edge of a read,
  2. encompasses more than one variant, or
  3. has the ImproperPair flag set 

The quality is incremented as follows:

distanceFromReadEdge = minimum distance from either end of the complete read context to the edge of the read  
baseQuality (SNV) = BASEQ at variant location  
baseQuality (MNV/Indel) = min BASEQ over core read context  
modifiedBaseQuality = min(baseQuality - `baseQualityFixedPenalty`, distanceFromReadEdge - `distanceFromReadEdgeFixedPenalty`) 

improperPairPenalty = `mapQualityImproperPaidPenalty` if improper pair flag set else 0  
distanceFromReference = number of somatic alterations to get to reference from the complete read context  
distanceFromReferencePenalty =  (distanceFromReference - 1) * `mapQualityAdditionalDistanceFromRefPenalty`  
modifiedMapQuality = MAPQ - `mapQualityFixedPenalty` - improperPairPenalty - distanceFromReferencePenalty  

matchQuality += max(0, min(modifiedMapQuality, modifiedBaseQuality))

If a `LENGTHENED` or `SHORTENED` jitter match is made we increment the jitter penalty as a function of the count of the repeat sequence in the microsatellite:

`JITTER_PENALTY` += `jitterPenalty` * max(0, repeatCount - `jitterMinRepeatCount`)

The final quality score also takes into account jitter and is calculated as:

`QUAL` =  matchQuality - `JITTER_PENALTY`

### Output

The outputs of this stage are found in the VCF as:

Field | Formula
---|---
RC_CNT\[0,1,2,3,4,5\] | Read Context Count \[`FULL`, `PARTIAL`, `CORE`, `REALIGNED`, `REFERENCE`, `TOTAL`\]
RC_JIT\[0,1,2\] | Read Context Jitter \[`SHORTENED`, `LENGTHENED`, `JITTER_PENALTY`\]
AD\[0,1\] | Allelic Depth \[`REFERENCE`, `FULL` + `PARTIAL` + `CORE` + `REALIGNED`\]
DP | Allelic Depth (=RC_CNT\[5\])
AF | Allelic Frequency AD\[1\] / DP
QUAL | Variant Quality

## 3. Hard Filters

To reduce processing time there are two hard filters that are eagerly applied at this stage. 

Filter | Default Value | Field
---|---|---
hard_min_tumor_qual | 1| `QUAL`
hard_min_tumor_alt_support |2| Normal `AD[1]`
 
These variants are excluded from this point onwards and have no further processing applied to them.
 
## 4. Normal Counts and Quality
Evidence in the normal is collected in same manner as step 2.

## 5. Soft Filters
TODO: CLEAN UP


Given evidence of the variants in the tumor and normal we apply somatic filters. 
The key principles behind the filters are ensuring sufficient support for the variant (minimum VAF and score) in the tumor sample 
and validating that the variant is highly unlikely to be present in the normal sample.
The filters are tiered to maximise sensitivity in regions of high prior likelihood for variants.
A hotspot panel of 10,000 specific variants are set to the highest sensitivity followed by medium sensitivity for a panel of cancer related 
gene exons and splice regions and more aggressive filtering genome wide to ensure a low false positive rate.   
The default filtering settings are:

The following filters are applied after collecting evidence about a variant. 
The filters are applied according to the `TIER` of the variant. 
The `TIER` can be one of `HOTSPOT`, `PANEL` or `WIDE` as determined by the supplied hotspot and panel locations. 

Filter  | Hotspot | Panel | Wide | Field
---|---|---|---|---
min_tumor_qual|35|100|150|`QUAL`
min_tumor_vaf|0.5%|1.5%|2.5%|`AF`
min_germline_depth|0|0|10 | ?
min_germline_depth_allosome|0|0|6 | ?
max_germline_vaf|10%|4%|4% | ?
max_germline_rel_base_qual|100%|4%|4% | ?


## 6. Phasing

Somatic variants can be phased using the complete read context with nearby germline variants or other somatic variants.

Phasing is interesting for somatic calling from 2 perspectives: 
  - understanding the somatic mutational mechanism which has led to the variant; and 
  - understanding the functional impact of the variation.
  
Regarding mechanism, multiple somatic cis-phased variants can frequently occur together with prominent mechanisms being 2 base MNVs (eg. CC>TT in UV Signatures and CC>AA in Lung smoking signatures) and micro-inversions (which are normally called as offsetting INS and DEL).

Phasing somatic and germline variants together can also aid in understanding somatic mechanisms - for example if the germline has a 6 base deletion in a microsatellite and the tumor has a 7 base deletion, then the likely somatic mechanism is a 1 base deletion. 

Phasing is also important for function impact particularly in 2 prominent cases: 
  - 2 nearby phased frameshift variants can lead to an inframe (and potentially activating) INDEL; and 
  - 2 phased synonymous SNVs in the same codon could potentially cause a nonsense or missense effect together.    

Two variants are considered phased if their read contexts are identical after adjusting for their relative position.
This is demonstrated in the example below where two SNVs share an identical sequence of bases.

<pre>
REF: CAACAATCGAACGATATAAATCTGAAA
A>T: CAACAATCGA<b>T</b>CGATACAATC
T>C:       TCGATCGATA<b>C</b>AAATCTGAAA
</pre>

Similarly, SNVs and INDELs may be phased together.

Any variants that are phased together will be given a shared local phase set (`LPS`) identifier.

Phasing variants opens up a number of algorithmic possibilities including MNV detection and de-duplication as explained below.


### INDEL De-duplication

While the read context is designed to capture a unique sequence of bases, it it sometimes possible that repeat sequences in the flanks of the read context coupled with an aligners alternate view on the same event can cause duplicate INDELs. 
If SAGE finds two phased INDELs of the same type at the same position where one is a subset of the other, then the longer is filtered with `dedup`.
 

## 7. MNV Handling

Phased SNVs separated by no more than one base may indicate the existence of a MNV. 
To confirm the existence of a MNV we re-examine the bam files, but only if:
  1. the SNVs are unfiltered;
  2. the resultant MNV is a hotspot; or
  3. one of the SNVs is unfiltered and the other is a germline variant (ie, soft-filtered by the germline criteria only).

Once the evidence is collected, the MNV is included if it is either unfiltered with <b> NO SUPPORT IN THE NORMAL </b> or it's a hotpot. 
In either case, the constituent SNVs are filtered with `merge`.   

If the MNV contains a germline variant (case 3 above) it is filtered with `germline_mnv` before being written to file and the constituent SNVs remain unchanged.
The purpose of this is to capture the impact of the germline variant together with the somatic variant. 

The following example shows unfiltered, phased, SNVs that produce the MNV `TCGA > CGGT`:
<pre>
REF: CAACAATCGATCGATATAAATCTGA
T>C: CAACAA<b>C</b>GGTTCGATATAAATC
C>G:  AACAAC<b>G</b>GTTCGATATAAATCT
A>T:    CAACGG<b>T</b>TCGATATAAATCTGA

MNV: CAACAA<b>CG</b>G<b>T</b>TCGATATAAATCTGA
</pre>


## 8. Output

There are two more 'hard' filters that are lazily applied at the end of the process just before writing to file. 
These do not save any processing time but do reduce the output file size. 

Filter | Default Value | Field
---|---|---
hard_min_tumor_qual_vcf | 30 | `QUAL`
hard_max_normal_alt_support |2| Normal `AD[1]`

Including the `hard_filter` flag will turn all the soft filters described above into (lazily applied) hard filters.
