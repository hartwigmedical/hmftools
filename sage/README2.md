
# Somatic Alterations in Genome (SAGE)
SAGE is a somatic SNV, MNV and small INDEL caller.


# Algorithm

SAGE parses the tumor bam twice. The first time to get a list of candidate variants and determine a unique read context of each variant. 
The second pass counts the number of full, partial and realigned matches of the read context.

Finally, SAGE parses the normal bam once to get germline statistics on the variants.  

## Variant Support

In the first pass of the tumor BAM, SAGE uses the `I` and `D` flag in the CIGAR to find INDELs. 
To find SNVs, SAGE compares the bases in every aligned region (flags `M`, `X` or `=`) with the provided reference genome. 
MNVs are not included at this stage, but are determined after phasing. Only records meeting the `min_map_quality` requirement are processed. 
There is no base quality requirement. 

The output of this stage is a set of candidate variants with counts of the ref support, alt support and read depth. 
This corresponds to the AD and DP fields found in the final VCF. 
Note these values do not contribute to the quality or VAF calculations. These are calculated in the next step. 


## Read Context

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
PARTIAL  | Core and at least one flank match fully at same position. Remaining flank matches but is truncated.
REALIGNED  | Core and both flanks match exactly but at a different position. 
SHORTENED | Core and both flanks match with the removal of one repeat. Can be at same or realigned position.
LENGTHENED | Core and both flanks match with the addition of one repeat. Can be at same or realigned position.


## Quality Score

Full and partial read context matches both contribute to the quality score. Realigned matches do not. Lengthened and shortened realigned matches subtract from the score.



## Phased Variants
Two variants are considered phased if their read contexts are identical after adjusting for their relative position.
This is demonstrated in the example below where two SNVs share an identical sequence of bases.

<pre>
REF: CAACAATCGAACGATATAAATCTGAAA
A>T: CAACAATCGA<b>T</b>CGATAAAATC
T>C:       TCGATCGATA<b>C</b>AAATCTGAAA
</pre>

Similarly, SNVs and INDELs may be phased together.

Any variants that are phased together will be given a shared local phase set (`LPS`) identifier.

Phasing variants opens up a number of algorithmic possibilities including MNV detection and de-duplication as explained below.

### MNV Detection

If two SNVs are separated by no more than one base, they will be merged together to form a MNV.  
For example, the following phased SNVs will me merged to give the resultant MNV `TCGA > CGGT`:
<pre>
REF: CAACAATCGATCGATATAAATCTGA
T>C: CAACAA<b>C</b>GGTTCGATATAAATC
C>G:  AACAAC<b>G</b>GTTCGATATAAATCT
A>T:    CAACGG<b>T</b>TCGATATAAATCTGA

MNV: CAACAA<b>CG</b>G<b>T</b>TCGATATAAATCTGA
</pre>

The smallest value of each of the SNV properties (ie ref support, alt support, read context full match etc) is used to construct the MNV.
The merged SNVs are filtered with `merge`.

### SNV De-duplication

By convention, INDELs start one base to the left of the actual insert or delete. 
If this base is also a SNV, this will result in 2 variants, 1 for the SNV and 1 for the INDEL, ie:

```
C > T
C > TCAA
```

If the two variants are phased, the SNV is superfluous, and is thus filtered with `dedup`. 


### INDEL De-duplication

While the read context is designed to capture a unique sequence of bases, it it sometimes possible that repeat sequences in the flanks of the read context coupled with an aligners alternate view on the same event can cause duplicate INDELs. 
If SAGE finds two phased INDELs of the same type at the same position where one is a subset of the other, then the longer is filtered with `dedup`.
 

## Output


