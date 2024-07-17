# Virus Interpreter

Virus Interpreter is a [VIRUSBreakend](https://pubmed.ncbi.nlm.nih.gov/33973999) post-process algo that takes in the final VIRUSBreakend
summary and adds annotation and interpretation and performs filtering for reporting. The algo writes a "virus.annotated.tsv" where every line is an
annotated viral presence from the VIRUSBreakend summary file.

### Annotation

Virus Interpreter picks the reference taxid that should be displayed in a report and performs a look-up in the taxonomy db to find the matching virus name.

### Interpretation

Virus Interpreter allows the mapping of the species taxid as defined in the reporting db ("HPV", "EBV" ,"MCV", "HBV" or "HHV-8"). 
Within the Hartwig pipeline, this configuration is used to map all clinically relevant HPV species to "HPV", 
which in turn is used to label the sample as "HPV positive" or "HPV negative".

### Reporting

Every virus detected by VIRUSBreakend is evaluated for reporting. For a virus to be reported, the following conditions must be met:
 - The VIRUSBreakend QC status must not be `LOW_VIRAL_COVERAGE`
 - The virus species should be potentially reportable (i.e. present in the reporting db)
 - Detection of at least 1 integration site into the tumor DNA
   - For "EBV", the following additional conditions must also be met:
     - percentage covered of the viral genome is greater than 90%
     - coverage of the virus is higher than expected clonal mean coverage
     - The QC status of sample should not be `FAIL_CONTAMINATION` or `FAIL_NO_TUMOR`
 - If VIRUSBreakend could not detect integration sites into the tumor DNA, the virus is still reportable if the following conditions are met:
   - percentage covered of the viral genome is greater than 90% 
   - coverage of the virus is higher than expected clonal mean coverage 
   - The QC status of sample should not be `FAIL_CONTAMINATION` or `FAIL_NO_TUMOR`

### Driver likelihood

Viruses that are potentially reportable (i.e. present in the reporting db) are annotated with a driver 
likelihood (HIGH/LOW) in the reporting db. Also, when a potentially reportable virus doesn't met the conditions of reportable the 
driver likelihood is set to UNKNOWN. 
All other viruses (i.e. the viruses that are not present in reporting db), will be annotated with driver likelihood UNKNOWN.

### Output data

Virus Interpreter produces a tsv file where every line (record) is an entry from the VIRUSBreakend summary file. 
The following fields are stored for each detected virus:

| Field                  | Description                                                                                                                                               |
|------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------|
| taxid                  | The reference taxid of the virus that is called with VIRUSBreakend                                                                                        |
| name                   | The name of the virus, matching with the taxid                                                                                                            |
| qcStatus               | The QC status as reported by VIRUSBreakend                                                                                                                |
| integrations           | The number of detected integrations of this virus into the sample genome as reported by VIRUSBreakend                                                     |
| interpretation         | The output of the interpretation step of Virus Interpreter                                                                                                |
| percentageCovered      | The percentage of the viral reference sequence that has been covered in the tumor sample as reported by VIRUSBreakend ("coverage" field of VIRUSBreakend) |
| meanCoverage           | The average coverage of the viral genome as reported by VIRUSBreakend  ("meanDepth" field of VIRUSBreakend)                                               |
| expectedClonalCoverage | The expected coverage assuming the virus is clonally integrated once in the tumor DNA                                                                     |
| reported               | A boolean indicating whether the detected viral presence is considered a driver                                                                           |
| driverLikelihood       | The driver likelihood of the virus as annotated in the reporting db                                                                                       |

## Version History and Download Links
 - Upcoming:
   - Add blacklist option 
 - [1.2](https://github.com/hartwigmedical/hmftools/releases/tag/virus-interpreter-v1.2)
   - Annotate the viruses with a driver likelihood
   - Version built on java11 JDK
 - [1.1](https://github.com/hartwigmedical/hmftools/releases/tag/virus-interpreter-v1.1)
   - New reporting strategy of viruses to report only clinical relevant viruses
 - [1.0](https://github.com/hartwigmedical/hmftools/releases/tag/virus-interpreter-v1.0)
   - Initial release. 