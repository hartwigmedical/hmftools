# Virus Interpreter

Virus Interpreter is a [VIRUSBreakend](https://pubmed.ncbi.nlm.nih.gov/33973999) post-process algo that takes in the final VIRUSBreakend
summary and adds annotation and interpretation and performs filtering for reporting. The algo writes a "virus.annotated.tsv" where every line is an
annotated viral presence from the VIRUSBreakend summary file.

### Annotation

Virus Interpreter picks the reference taxid that should be displayed in a report and performs a look-up in the taxonomy db to find the matching virus name.

### Interpretation

Virus Interpreter allows the mapping of the species taxid which are definied into the reporting db ("HPV", "EBV" ,"MCV", "HBV" or "HHV-8"). 
Within the Hartwig pipeline this configuration is used to map all clinically relevant HPV species to "HPV" 
which in turn is used to label patients as "HPV positive" or "HPV negative".

### Reporting

Every virus found by VIRUSBreakend is evaluated for reporting. For a virus to be reported, the following conditions need to be met:
 - QC status of sample should not contains `FAIL`
 - The VIRUSBreakend QC status must not be `LOW_VIRAL_COVERAGE`
 - The virus should be present in the reporting db with the conditions for reporting
 - VIRUSBreakend must have found at least 1 integration site into the tumor DNA for "HPV", "MCV", "HBV" or "HHV-8"
   - For "EBV" next to the at least 1 integration site the following conditions should extend with: 
     - percentage covered of the virus should be greater than 90%
     - coverage of virus should be greater than the expected clonal coverage
 - VIRUSBreakend has none integration sites into the tumor DNA for "HPV", "MCV", "HBV", "EBV" or "HHV-8" and the conditions should extend with: 
   - percentage covered of the virus should be greater than 90% 
   - coverage of virus should be greater than the expected clonal coverage 
   
### Output data

Virus Interpreter produces a tsv file where every line (record) is an entry from the VIRUSBreakend summary file. 
The following fields are stored per viral presence:

Field | Description 
---|---
taxid | The reference taxid of the virus that is called with VIRUSBreakend
name | The name of the virus, matching with the taxid
qcStatus | The QC status as reported by VIRUSBreakend
integrations | The number of detected integrations of this virus into the sample genome as reported by VIRUSBreakend
interpretation | The output of the interpretation step of Virus Interpreter
percentageCovered | The percentage of the viral reference sequence that has been covered in the tumor sample as reported by VIRUSBreakend
coverage | The mean coverage of the virus as reported by VIRUSBreakend  //TODO: improve
expectedClonalMeanCoverage | The expected coverage assuming the virus is clonally integrated once in the tumor DNA 
reported | A boolean indicating whether the detected viral presence is considered a driver

 ## Version History and Download Links
 - [1.1] (coming)
   - New reporting strategy of viruses to report only clinical relevant viruses
 - [1.0](https://github.com/hartwigmedical/hmftools/releases/tag/virus-interpreter-v1.0)
   - Initial release. 