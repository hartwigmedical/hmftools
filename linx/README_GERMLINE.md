
## Germline Linx
Linx may be run in a germline mode which chains on proximity only and does not rely on copy number.  It is optimised to detect germline pathogenic structural variants.
## Example Usage
This is a typical command to run LINX for a single sample from PURPLE output.

```
java -jar linx.jar 
TO DO!!!!!!!!!!!
```
## Algorithm
For the Germline, LINX is run in a limited mode  where we perform only the following operations
- Artefact filtering
- Proximity clustering
- LINE annotation
- Chaining of assembled links only
- Breakend annotation of gene impact for genes in the driver gene panel.  Note that PON filtered SV are also annotated with gene impacts for germline SV.    

The same logic is used to annotate breakends as disruptive as for somatic SVs, with special additional logic that a deletion (or single breakend with a candidate mapping that could indicate a deletion) that completely encompasses a gene is also considered disruptive (this is not handled in the somatic logic as it is treated as a copy number event).  

LINX additionally filters any deletion which exactly matches exon boundaries on both ends as ‘PSEUDOGENE’

Disruptive breakends are marked as ‘reportable’ and added to the driver catalog if:
- Filter = “PASS”
- The variant is biologically plausible in the germline.   ie. Resolved Type = “DEL”, “DUP”,  “RECIPROCAL_INVERSION” or “RECIPROCAL_TRANSLOCATION”, including synthetic events.   We also report any variant in a more complex cluster that is of type “DEL” or “DUP” .  For “DEL” and “DUP” we also require that the length of the event is < 3M bases.   “SGL” breakends are also allowed if a candidate alignment exists that would make it a reportable DEL or DUP <500Kb bases.  
The impacted gene has been configured as report_germline_disruptions=”TRUE” in the driver gene panel.
- CohortFrequency (ie PonCount) < 3

## Output

Field | Description
--|--
SampleId|Id of sample
ChromosomeStart|As reported by GRIDSS
PositionEnd|As reported by GRIDSS
OrientationStart|As reported by GRIDSS
ChromosomeEnd|As reported by GRIDSS
PositionStart|As reported by GRIDSS
OrientationEnd|As reported by GRIDSS
Gene|Genes impacted
Type|‘INS’,’DEL’,’DUP’,’SGL’,’INV’ or ‘BND’
Filter|‘PASS’ or ‘PON’
Event|Gridss identifier
QualScore|Quality score for reference sample only
GermlineVariantFragments|Fragments supporting the germline SV
GermlineReferenceFragmentsStart|Fragments supporting the reference at the start position
GermlineReferenceFragmentsEnd|Fragments supporting the reference at the end position
TumorVariantFragments|Fragments supporting the somatic SV
TumorReferenceFragmentsStart|Fragments supporting the somatic SV at the start position
TumorReferenceFragmentsEnd|Fragments supporting the somatic SV at the end position
ResolvedType|Cluster classification
LinkedByStart|Id of any linked assemblies or linked double stranded breaks
LinkedByEnd|Id of any linked assemblies or linked double stranded breaks
insertSequenceAlignments|As reported by GRIDSS
insertSequenceRepeatClass|As reported by GRIDSS
insertSequenceRepeatType|As reported by GRIDSS
CohortFrequency|PON Count (if variant is found in PON)
Reported|T/F



