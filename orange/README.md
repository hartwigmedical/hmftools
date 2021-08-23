# Oncogenic Results of Analyzing the Genome

ORANGE summarizes the key outputs from all algorithms in the Hartwig suite into a single PDF:
 1. The algo depends exclusively on config and data produced by the [Hartwig platinum pipeline](https://github.com/hartwigmedical/platinum) 
 and hence can always be run as final step without any additional local data or config required. 
 1. The algo intends to combine RNA and DNA data to present an integrated DNA/RNA analysis of a tumor sample. 
 1. Everything that is labeled as a driver by any of the Hartwig algo's is displayed in the PDF along with the driver likelihood. 
 This effectively means that everything reported by [patient-reporter](../patient-reporter) is present in the ORANGE pdf.
 1. An additional exhaustive WGS scan is performed for anything interesting that may be potentially relevant but not picked up as a driver.
 Details of what is considered interesting are described in below.
 1. A broad range of QC measures and plots is displayed which should make it clear what the quality of the tumor sample is without having to 
 dig into additional outputs of QC steps.     
 
An example report based on the publicly available melanoma cell line COLO829 can be found [here](src/main/resources/Test.orange.pdf).

Note that neither this readme nor the report itself contains any documentation about the Hartwig algorithms and output. For questions in 
this area please refer to the specific algorithm documentation present on [https://github.com/hartwigmedical/hmftools](https://github.com/hartwigmedical/hmftools)

The front page of the ORANGE report lists all high-level stats about the sample along with genome-wide visualisations of all mutations and 
SNV/Indel clonality. In addition to this front page, the following chapters are generated in the ORANGE report:
 
  - [Clinical Evidence](#clinical-evidence): What genomic evidence has been found in favor of, or against, specific treatments?
  - [Somatic Findings](#somatic-findings): What potentially relevant mutations have been found in the tumor specifically?
  - [Germline Findings](#germline-findings): What potentially relevant mutations have been found in the germline DNA? 
  - [Immunology](#immunology): What can we tell about the immunogenicity of the tumor sample?
  - [Cohort Comparison](#cohort-comparison): How do the various properties of this tumor compare to existing cancer cohorts?
  - [Quality Control](#quality-control): Various stats and graphs regarding the quality of the data and interpretation thereof. 

### Optional Configuration

Argument | Description
---|---
disable_germline  | If set, disables the germline findings chapter and removes clinical evidence on germline events.
max_evidence_level | If set, filters evidence down to this level. For example, if "B" is passed as a parameter, only treatments with at least A or B level evidence are displayed in the clinical evidence chapter of the report. Do note that the front page always lists the count of all evidence present, regardless of this filter setting.

### Clinical Evidence 

The following algo is used to render clinical evidence in the ORANGE report based on [PROTECT](../protect) output:
 1. All evidence is considered for display, including evidence based on non-drivers or evidence that is normally filtered for reporting by PROTECT.
 1. Evidence is split between trials and non-trials which are further split up based on on/off label. 
 1. Evidence is grouped by treatment and split up between responsive and resistance evidence.
 1. Evidence is filtered for optional configuration (germline, max_reporting_level)
 
Effectively this means that the evidence is _exhaustive_ but not necessarily all _applicable_. 

### Somatic Findings 

In addition to all somatic drivers (SNVs/Indels, copy numbers, structural variants and fusions) the following is considered potentially
interesting and added to the report:
 1. Other potentially relevant variants
    * Variants that are hotspots but not part of the reporting gene panel.
    * Variants which have clinical evidence but are not part of the reporting gene panel.
    * Coding variants that are not reported but are phased with variants that are reported.
    * Variants that are considered relevant for tumor type classification according to CUPPA.
 1. Other regions with amps or autosomal losses:
    * Any chromosomal band location with amplification or loss is picked for reporting
    * For a band with at least one gene amplified, the gene with the highest copy number is picked.
    * For a band with a loss that has no losses reported in this band a random gene is picked
    * A maximum of 10 additional gains and 10 additional losses are reported as potentially interesting. 
 1. Other potentially relevant fusions:
    * Any fusion that is not reported and has a reported type other than NONE is picked. 
    * Any fusion with clinical evidence is picked. 
    * A maximum of 10 additional fusions is reported as potentially interesting.
 1. Other viral presence
    * Any viral presence that is not otherwise reported is reported through this route. 
    
### Germline Findings

In addition to all germline SNV/Indel tumor drivers determined by [PURPLE](../purple), the following is added to the report:
 1. Other potentially relevant variants
    * Any hotspots that are not configured to be reported.
    * Any hotspots that are filtered based on quality.
    
The germline CN aberrations are determined by [PURPLE](../purple) and include aberrations such as klinefelter or trisomy X. 

### Immunology

The immunology chapter is work-in-progress and will report on various immunology properties of the tumor sample.

### Cohort Comparison

The cohort comparison reports all the properties of a tumor sample that [CUPPA](../cuppa) considers for determining tumor type. The cohort
comparison displays the prevalence of the tumor's properties with respect to the cohorts that CUPPA could potentially assign the sample to:
 1. Genomic position distribution of SNVs and their tri-nucleotide signature
 1. Sample traits of the tumor (for example, number of LINE insertions)
 1. (Driver) features of the tumor.
 
Do note that RNA features and cohort comparison thereof are only included if platinum was run in combined DNA/RNA mode.  
 
### Quality Control

While the overall QC is usually reliable, the quality control chapter displays all stats and plots relevant for interpreting the overall QC:
 1. The high-level QC from [PURPLE](../purple)
 1. Various details from the sample's flagstats and coverage stats
 1. Various plots from [PURPLE](../purple)  
 1. BQR plots from both reference and tumor sample from [SAGE](../sage)


