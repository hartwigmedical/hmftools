
# Somatic Alterations in Genome (SAGE)
SAGE is a somatic SNV, MNV and small INDEL caller.


## Distinct Read Context

The read context is the distinct set of bases surrounding the variant after accounting for any repeats and microhomology in the read. 
In this context, a repeat is defined as having 1 - 10 bases repeated at least 2 times. 

For an SNV in a non-repeat sequence this will just be the single alternate base, 'A'. 
For an SNV in a repeat, the entire repeat will be included as well as one base on either side, eg 'TAAAAC'.

A DEL will always include the bases on either side of the deleted sequence. 
If the deleted read sequence is part of a microhomology or repeat sequence, this will also be included in the read context.


The importance of capturing the microhomology in the read context is demonstrated below. 
In the following example there is a delete of 4 bases that will nominally be left aligned at position 7 as AAAAC > A. 
     

Given a DEL of AAAAC > A

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
ALT:   GTCT<b>CAAAAACAA    ACAAACAAT</b>AAAAAAC
ALT:   GTCT<b>CAAAAACAAA    CAAACAAT</b>AAAAAAC
ALT:   GTCT<b>CAAAAACAAAC    AAACAAT</b>AAAAAAC
ALT:   GTCT<b>CAAAAACAAACA    AACAAT</b>AAAAAAC
ALT:   GTCT<b>CAAAAACAAACAA    ACAAT</b>AAAAAAC
ALT:   GTCT<b>CAAAAACAAACAAA    CAAT</b>AAAAAAC
ALT:   GTCT<b>CAAAAACAAACAAAC    AAT</b>AAAAAAC
ALT:   GTCT<b>CAAAAACAAACAAACA    AT</b>AAAAAAC
ALT:   GTCT<b>CAAAAACAAACAAACAA    T</b>AAAAAAC

REF:   GTCTCAAAAACAAACAAACAAACAATAAAAAAC 
ALT:   GTCTCAAA    AACAAACAAACAATAAAAAAC
ALT:   GTCTCAAAA    ACAAACAAACAATAAAAAAC
ALT:   GTCTCAAAAA    CAAACAAACAATAAAAAAC
ALT:   GTCTCAAAAAC    AAACAAACAATAAAAAAC
ALT:   GTCTCAAAAACA    AACAAACAATAAAAAAC
ALT:   GTCTCAAAAACAA    ACAAACAATAAAAAAC
ALT:   GTCTCAAAAACAAA    CAAACAATAAAAAAC
ALT:   GTCTCAAAAACAAAC    AAACAATAAAAAAC
ALT:   GTCTCAAAAACAAACA    AACAATAAAAAAC
ALT:   GTCTCAAAAACAAACAA    ACAATAAAAAAC
ALT:   GTCTCAAAAACAAACAAA    CAATAAAAAAC
ALT:   GTCTCAAAAACAAACAAAC    AATAAAAAAC
ALT:   GTCTCAAAAACAAACAAACA    ATAAAAAAC
ALT:   GTCTCAAAAACAAACAAACAA    TAAAAAAC
</pre>






Ref | Alt | ReadContext | Repeat | Microhomology
------------ | -------------| -------------| -------------| -------------
 C | A | C | NA | NA
 C | A | TAAAAAC | Ax4 | NA
 CCGT | C | CT | NA | NA
 AAAAC | A | CAAAAACAAACAAACAAT | Ax5 | AAAC
 