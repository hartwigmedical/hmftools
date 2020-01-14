
# Somatic Alterations in Genome (SAGE)
SAGE is a somatic SNV, MNV and small INDEL caller.


# Algorithm

SAGE parses the tumor bam twice. The first to get a list of candidate variants and determine a unique read context of each variant. 
The second pass counts the number of full, partial and realigned matches of the read context.

Finally, SAGE parses the normal bam once to get germline statistics on the variants.  

## Tumor First Pass - Candidates

In the first pass of the tumor BAM, SAGE uses the `I` and `D` flag in the CIGAR to find INDELs and compares the bases in every aligned 
region (flags `M`, `X` or `=`) with the provided reference genome to find SNVs. 
Note that there are no base quality or mapping quality requirements when looking for candidates.

SAGE tallies the raw ref/alt support and base quality and collects the read contexts of each variant.
Once finished, each variant is assigned its most frequently found read context as its primary one. 
If a variant does not have at least one complete read context (including flanks) it is discarded.
All remaining variants are then considered candidates for processing in the second pass. 

The variants at this stage have the following properties available in the VCF:

Field | Description
---|---
RC | Read context (core only without flanks)
RDP | Raw depth
RAD\[0,1\] | Raw allelic depth \[Ref,Alt\]
RABQ\[0,1\] | Raw allelic base quality \[Ref,Alt\]

Note that these values do NOT contribute to the AD, DP, QUAL or AF fields. These are calculated in the second pass. 

## Tumor Second Pass - Counts and Quality


# Read Context

The read context is a sequence of bases comprised of a core flanked on either side by an additional 25 bases. 
A read context must be complete to be eligible as the primary read context of a variant. 
If a single flank is incomplete (such as when a variant is too close to a read edge), the read context can partially match with a complete read context.
If both flanks are incomplete (such as when the majority of a read is in a repeat sequence), the read context cannot match with another read context.


### Read Context Core

The read context core is the distinct set of bases surrounding the variant after accounting for any repeats and microhomology in the read sequence (not ref sequence). 
In this context, a repeat is defined as having 1 - 10 bases repeated at least 2 times. 

For a SNV in a non-repeat sequence this will just be the single alternate base. 
For a SNV in a repeat, the entire repeat will be included as well as one base on either side, eg 'TAAAAC'.

A DEL will always include the bases on either side of the deleted sequence. 
If the deleted read sequence is part of a microhomology or repeat sequence, this will also be included in the read context.

An INSERT will always include the base to the left of the insert as well as the new sequence. 
As with a DEL, the read context will be extended to include any repeats and/or microhomology.

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

### Read Context Counts

Variant quality and VAF are determined using counts of the read context rather than the variant support.  

Type  | Description
---|---
FULL | Core and both flanks match all match at same position. 
PARTIAL | Core and at least one flank match fully at same position. Remaining flank matches but is truncated.
CORE | Core matches but neither flank does.
REALIGNED | Core and both flanks match exactly but at offset from the expected position. 
REFERENCE | Read matches the reference at the position. 
TOTAL | Total number of reads that cover the read context (excluding flanks).


These values are found in the `RCC` field in the VCF output along with the following derived fields:

Field | Formula
---|---
`RC_CNT\[0,1,2,3,4,5\]` | \[FULL, PARTIAL, CORE, REALIGNED, REFERENCE, TOTAL\]
`AD\[0,1\]` | \[REFERENCE, FULL + PARTIAL + CORE + REALIGNED\]
`DP` | TOTAL
`AF` | AD\[1\] / DP


SHORTENED | Core and both flanks match with the removal of one repeat. Can be at same or realigned position.
LENGTHENED | Core and both flanks match with the addition of one repeat. Can be at same or realigned position.


## Quality Score

Full and partial read context matches both contribute to the quality score. Realigned matches do not. Lengthened and shortened realigned matches subtract from the score.



## Phased Variants
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

### MNV Detection

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


### INDEL De-duplication

While the read context is designed to capture a unique sequence of bases, it it sometimes possible that repeat sequences in the flanks of the read context coupled with an aligners alternate view on the same event can cause duplicate INDELs. 
If SAGE finds two phased INDELs of the same type at the same position where one is a subset of the other, then the longer is filtered with `dedup`.
 

## Output


# Filters

## Hard Filters

To reduce processing time there are two hard filters that are eagerly applied.  

Filter | Default Value | Field
---|---|---
hard_min_tumor_qual | 1| `QUAL`
hard_min_tumor_alt_support |2| Normal `AD[1]`

There are two more 'hard' filters that are lazily applied at the end of the process just before writing to file. 
These do not save any processing time but do reduce the output file size. 

Filter | Default Value | Field
---|---|---
hard_min_tumor_qual_vcf | 30 | `QUAL`
hard_max_normal_alt_support |2| Normal `AD[1]`

Including the `hard_filter` flag will turn all the soft filters below into (lazily applied) hard filters.

## Soft Filters

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



