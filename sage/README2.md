
# Somatic Alterations in Genome (SAGE)
SAGE is a somatic SNV, MNV and small INDEL caller.

## Candidate Detection





## Distinct Read Context

The read context is the distinct set of bases surrounding the variant after accounting for any repeats and microhomology in the read sequence (not ref sequence). 
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

TODO: explain sides/wings

## Quality Score

Full and partial read context matches both contribute to the quality score. Realigned matches do not. Lengthened and shortened realigned matches subtract from the score.



## Phased Variants
Two variants are considered phased if their read contexts are identical after adjusting for their relative position.
This is demonstrated in the example below where two SNVs shared an identical sequence of bases.

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
REF: CAACAATCGATCGATATAAATCTGAA
T>C: CAACAA<b>C</b>GGTTCGATATAAATC
C>G:  AACAAC<b>G</b>GTTCGATATAAATCT
A>T:    CAACGG<b>T</b>TCGATATAAATCTGA
MNV: CAACAA<b>CG</b>G<b>T</b>TCGATATAAATC
</pre>

The smallest value of each of the SNV properties (ie ref support, alt support, read context full match etc) is used to construct the MNV.
The merged SNVs are filtered with a flag of 'mnv'.

## Output


