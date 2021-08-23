# Oncogenic Results of Analyzing the Genome

ORANGE summarizes the key outputs from all algorithms in the Hartwig suite into a single PDF:
 1. The algo depends exclusively on config and data produced by the [Hartwig platinum pipeline](https://github.com/hartwigmedical/platinum) 
 and hence can always be run as final step without any local data or config required. 
 1. The algo intends to combine RNA and DNA data to present an integrated DNA/RNA analysis of a tumor sample. 
 1. Everything that is labeled as a driver by any of the Hartwig algo's is displayed in the PDF along with the driver likelihood. 
 This effectively means that everything reported by [patient-reporter](../patient-reporter) is present in the ORANGE pdf.
 1. An additional exhaustive WGS scan is performed for anything interesting that may be potentially relevant but not picked up as a driver.
 Details of what is considered interesting are described in below. 
 1. A broad range of QC measures and plots is added which should make it clear what the quality of the tumor sample is without having to 
 dig into additional outputs of QC steps.     
 
An example report based on the publicly available melanoma cell line COLO829 can be found [here](src/main/resources/Test.orange.pdf).

Note that neither this readme or the report itself contains any documentation about the Hartwig algorithms and output. For questions in 
this area please refer to the specific algorithm documentation present on [https://github.com/hartwigmedical/hmftools](https://github.com/hartwigmedical/hmftools)

The front page of the ORANGE report lists all high-level stats about the sample along with genome-wide visualisations of all mutations and 
SNV/Indel clonality. In addition to this front page, the following chapters are generated in the ORANGE report:
 
  - [Clinical Evidence](#clinical-evidence): What genomic evidence has been found in favor of, or against, specific treatments?
  - [Somatic Findings](#somatic-findings): What potentially relevant mutations have been found in the tumor specifically?
  - [Germline Findings](#germline-findings): What potentially relevant mutations have been found in the germline DNA? 
  - [Immunology](#immunology): What can we tell about the immunogenicity of the tumor sample?
  - [Cohort Comparison](#cohort-comparison): How does the various properties of this tumor compare to existing cancer cohorts?
  - [Quality Control](#quality-control): Various stats and graphs regarding the quality of the data and interpretation thereof. 

### Optional Configuration

Argument | Description
---|---
disable_germline  | If set, disables the germline findings chapter and removes clinical evidence on germline events.
max_evidence_level | If set, filters evidence down to this level. For example, if "B" is passed as a parameter, only drugs with at least A or B level evidence are displayed in the clinical evidence section of the report. Do note that the front page always lists the count of all evidence present, regardless of this filter setting.

### Clinical Evidence 

### Somatic Findings 

### Germline Findings

### Immunology

### Cohort Comparison

### Quality Control



