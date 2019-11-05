
# Somatic Alterations in Genome (SAGE)
SAGE is a somatic SNV, MNV and small INDEL caller.


## Distinct Read Context

The read context is the distinct set of bases surrounding the variant after accounting for any repeats and microhomology in the read. 
In this context, a repeat is defined as having 1 - 10 bases repeated at least 2 times. 

For an SNV in a non-repeat sequence this will just be the single alternate base, 'A'. 
For an SNV in a repeat, the entire repeat will be included as well as one base on either side, eg 'TAAAAC'.

A DEL will always include the bases on either side of the deleted sequence. 
If the deleted read sequence is part of a microhomology or repeat sequence, this will also be included in the read context.


Ref | Alt | ReadContext | Repeat | Microhomology
------------ | -------------| -------------| -------------| -------------
 C | A | C | NA | NA
 C | A | TAAAAAC | Ax4 | NA
 CCGT | C | CT | NA | NA
 AAAAC | A | CAAAAACAAACAAACAAT | Ax5 | AAAC
 