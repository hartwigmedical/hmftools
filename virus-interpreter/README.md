# Virus Interpreter

Virus Interpreter is a [VIRUSBreakend](https://pubmed.ncbi.nlm.nih.gov/33973999) post-process algo that takes in the final VIRUSBreakend
summary and adds annotation and interpretation and performs filtering for reporting. The algo writes a "virus.annotated.tsv" where every line is an
annotated viral presence from the VIRUSBreakend summary file.

### Annotation

Virus Interpreter picks the reference taxid that should be displayed in a report and performs a look-up in the taxonomy db to find the matching virus name.

### Interpretation

Virus Interpreter allows the mapping of any species taxid to either "HPV", "EBV" or "MCV". 
Within the Hartwig pipeline this configuration is used to map all clinically relevant HPV species to "HPV" 
which in turn is used to label patients as "HPV positive" or "HPV negative".

### Reporting

Every virus found by VIRUSBreakend is evaluated for reporting. For a virus to be reported, the following conditions need to be met:
 - VIRUSBreakend must have found at least 1 integration site into the tumor DNA
 - The VIRUSBreakend QC status must not be `LOW_VIRAL_COVERAGE`
 - The virus must not be blacklisted.
 
The blacklist is configurable and used in the Hartwig pipeline to filter any forms of HIV from getting reported.

### Output data

Virus interpreter produces a tsv file where every line (record) is an entry from the VIRUSBreakend summary file. 
The following fields are stored per viral presence:

Field | Description 
---|---
taxid | The taxid of the virus that is reported from VIRUSBreakend
name | The name of the virus, matching with the taxid
qcStatus | The QC status as produced by VIRUSBreakend
integrations | The number of integrations of this virus into the sample genome as reported by VIRUSBreakend
interpretation | The output of the Interpretation step of Virus Interpreter
reported | A boolean indicated this viral presence will be reported 


 ## Version History and Download Links
 - [1.0](https://github.com/hartwigmedical/hmftools/releases/tag/virus-interpreter-v1.0)
   - Initial release. 