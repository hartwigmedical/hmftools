
# Somatic Alterations in Genome (SAGE)
SAGE is a somatic SNV, MNV and small INDEL caller.


## Distinct Read Context

The read context is the distinct set of bases surrounding the variant after accounting for any repeats and microhomology in the read sequence (not the ref genome). 
In this context, a repeat is defined as having 1 - 10 bases repeated at least 2 times. 

For a SNV in a non-repeat sequence this will just be the single alternate base, 'A'. 
For a SNV in a repeat, the entire repeat will be included as well as one base on either side, eg 'TAAAAC'.

A DEL will always include the bases on either side of the deleted sequence. 
If the deleted read sequence is part of a microhomology or repeat sequence, this will also be included in the read context.

An INSERT will always include the base to the left of the insert as well as the new sequence. 
As with a DEL, the read context will be extended to include any repeats and/or microhomology.


The importance of capturing the microhomology is demonstrated in the following example. This delete of 4 bases is nominally left aligned at 
position 7 as AAAAC > A but can equally be  represented as 8:AAACA > A, 9:AACAA > A, 10: ACAAA > A, 11: CAAAC > C etc. 
Using a read context of `CAAAAACAAACAAACAAT` which spans the microhomology will match (shown bolded) every alt but not the ref:

The ref and read sequences appear as:

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






Ref | Alt | ReadContext | Repeat | Microhomology
------------ | -------------| -------------| -------------| -------------
 C | A | C | NA | NA
 C | A | TAAAAAC | Ax4 | NA
 CCGT | C | CT | NA | NA
 AAAAC | A | CAAAAACAAACAAACAAT | Ax5 | AAAC
 