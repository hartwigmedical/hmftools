# Cuppa

## Overview
Cuppa is a tool that weighs multiple genomic features observed from WGS & WTS to predict the tissue of origin of a tumor sample. 
It is intended to allow both prediction of primary tumor location for Cancer of Unknown Primary (CUP) samples and identify potential errors in clinical information.

The key inputs to Cuppa are:
- a collection of reference data files which define rates of various DNA and RNA features and characteristics
- the measurements of those same DNA and RNA features for the sample or samples being evaluated

Cuppa then runs 1 or more of the following 'classifiers' to make its assessment of a sample:

### DNA Classifiers
There are 3 DNA-based classifiers
- SNV trinucleotide counts - calculated using a pairwise test against each reference sample's SNV counts
- SNV genomic position frequencies - frequency of SNVs in 500K positional buckets, calculated using a test against per-cancer-type counts
- features - drawing on any other somatic event such as gene fusions, drivers, structural variant characteristics, viral insertions or 
other recurrent events, and using a likelihood-based calculation model

### RNA Classifiers
There are 2 RNA-based classifiers
- Gene expression - calculated using a pairwise test against each reference sample's gene expression TPMs
- Alternate splice junctions - using splice junctions undocumented by Ensembl, and then calculated from per-cancer-type average fragment support for each alternate splice junction

Each classifier is optional and run if enabled by config.

## Reference Data
Cuppa can be used to produce reference data files itself - see instructions below.

Otherwise they can be obtained from the HMF Resources page  [HMFTools-Resources](https://resources.hartwigmedicalfoundation.nl/).

Cuppa will attempt to load these files from a reference directory supplied in its command-line arguments. The files are:

Data Type | Source | Filename | Description 
---|---|---|---
Sample Set | ALL | cup_ref_sample_data | List of the reference samples - SampleId and CancerType
Features | DNA | cup_ref_feature_prev.csv | Prevalence of fusions, viral insertions, drivers and known INDELs per cancer type
Features | DNA cup_ref_driver_avg.csv | Average driver counts per cancer types 
Sample Traits | DNA | cup_ref_sample_trait_percentiles.csv | Percentiles for purity, ploidy and MS Indels per cancer type
Sample Traits | DNA | cup_ref_sample_trait_rates.csv | Whole genome duplicate and gender rates per cancer type
Sample Traits | DNA | cup_ref_gender_rates.csv | Optional overrides for expected rates of gender per cancer type 
SNVs | DNA | cup_ref_snv_counts.csv | Matrix of trinucleotide buckets (in rows) for each sample (columns)
SNVs | DNA | cup_ref_sample_pos_freq_counts.csv | Matrix of genomic position frequency counts (in rows) for each sample (columns) 
SNVs | DNA | cup_ref_sig_percentiles.csv | Percentiles of each signature per cancer type
SVs | DNA | cup_ref_sv_percentiles.csv | Percentiles of key structural variant characteristics per cancer type


## Data Sources
Sample data for both the reference data and the samples being evaluated can be sourced from any of 3 different inputs:
- a MySQL HMF Patients database
- flat-files from a HMF pipeline run (Purple, Linx, Isofox, Sage and GRIDSS)
- generic flat files

### Database

Data Type | Table 
---|---
Features | svFusion, viralInsertion, driverCatalog and somaticVariant (for known INDELs)
SNVs | somaticVariant
Traits | purity

### Pipeline Files

Data Type | File Details 
---|---
Features | Linx fusion, viral insert and driver catalog files, Purple somatic VCF
SNVs | Purple somatic VCF
Traits | Purple purity file

### Cohort Files
Cuppa also accepts CSV inputs files conforming to the same file format as produced when generating reference data as described below.
This can be handy and much more efficient when testing a large cohort in a single run.

Data Type | Cohort Filename | Fields & Comments
---|---|---
Features | cup_ref_cohort_feature_data.csv | SampleId,Name,Type(DRIVER,FUSION,VIRUS or INDEL),Likelihood,ExtraInfo
SNVs | cup_ref_snv_counts.csv | Matrix of trinucleotide counts in rows, SampleIds in columns 
SNVs | cup_ref_sample_pos_freq_counts.csv | Matrix of genomic position frequency counts in rows, SampleIds in columns
SNVs | cup_ref_cohort_signature_data.csv | SampleId,Signature,AllocationPercentage 
Traits | cup_ref_cohort_traits_data.csv | SampleId,Gender,WholeGenomeDuplication,Purity,Ploidy,MsIndelsPerMb,ChordHrd
SVs | cup_ref_cohort_sv_data.csv | SampleId,LINE,SIMPLE_DEL_20KB_1MB,SIMPLE_DUP_32B_200B,SIMPLE_DUP_100KB_5MB,MAX_COMPLEX_SIZE,TELOMERIC_SGL


## Running Cuppa

### Mandatory Arguments

Argument | Description 
---|---
sample_data | Sample ID
ref_data_dir | Reference data directory
sample_data_dir | Sample data directory containing Linx and Purple files
sample_sv_file | Sample structural variant VCF
sample_somatic_file | Sample somatic variant VCF
output_dir | Path to write sample Cuppa output

### Optional Arguments

Argument | Description 
---|---
categories | By default Cuppa will run all DNA categories. A subset can be specified as a ';' list from SNV, SV, SAMPLE_TRAIT, FEATURE
write_similarities | Write the top-20 cosine similarities for SNVs

### Example using pipeline files
```
java -jar cuppa_jar \
    -categories DNA \
    -ref_data_dir /reference_data_dir/ \
    -sample_data SAMPLE_ID \
    -sample_data_dir /sample_pipeline_files_dir/ \
    -sample_sv_file /sample_sv_vcf_file/ \
    -sample_somatic_vcf /sample_snv_vcf_file/ \
    -output_dir /output_dir/ \
```

### Example using database
```
java -jar cuppa_jar \
    -categories DNA \
    -ref_data_dir /reference_data_dir/ \
    -sample_data SAMPLE_ID \
    -db_url DB_URL -db_user DB_USER -db_pass DB_PASS \
    -output_dir /output_dir/ \
```

## Reference data generation
Cuppa can be used to generate the reference files used for subsequently evaluating samples.

To generate reference data, the key input is a list of sample IDs and their designated cancer-type classification. 
The unique set of cancer types defined in the reference data will then drive the evaluation of each sample subsequently tested by Cuppa.

### Mandatory Arguments
Argument | Description
---|---
ref_sample_data | a CSV file of SampleId,CancerType for each reference sample
db_url,db_user,db_pass | Connection details to MySQL HMF patients DB

### Optional Arguments
Argument | Description
---|---
gender_rates | Optional, expected rate of cancer for specific cancer types in form 'CancerType1;FemaleRate1;MaleRate1,CancerType2;FemaleRate2;MaleRate2'
feature_override_file | Override specific feature prevalences
write_cohort_files | Rewrites the reference sample's data for each category for subsequent optimised use in Cuppa 

An example command to generate reference data for a cohort is shown below. 
The file 'cup_ref_sample_data.csv' is a CSV file of SampleId,CancerType for each reference sample. 

```
java -cp cuppa.jar com.hartwig.hmftools.cup.ref.RefDataBuilder \ 
  -ref_sample_data_file cup_ref_sample_data.csv 
  -db_url DB_URL -db_user DB_USER -db_pass DB_PASSWORD \ 
  -gender_rates "Breast;1;0.1" \
  -write_cohort_files \
  -output_dir /output_dir \
```

If cohort files are already available, then they can be used to generated reference data instead of sourcing data from the database:
```
java -cp cuppa.jar com.hartwig.hmftools.cup.ref.RefDataBuilder \ 
  -ref_sample_data_file cup_ref_sample_data.csv 
  -ref_snv_counts_file cup_ref_snv_counts.csv \
  -ref_sample_snv_pos_freq_file cup_ref_sample_pos_freq_counts.csv \
  -ref_sample_traits_file cup_ref_cohort_traits_data.csv \
  -ref_sig_contribs_file cup_ref_cohort_signature_data.csv \
  -ref_sv_data_file cup_ref_cohort_sv_data.csv \
  -ref_features_file cup_ref_cohort_feature_data.csv \ 
  -output_dir /output_dir \
```

## Output Files
Cuppa writes a probability for each feature or characteristic per cancer type. From this is computes an overall 'COMBINED' probability per cancer type.

The output file has these fields

Field | Description
---|---
SampleId | Sample being evaluated
Category | SNV, SV, SAMPLE_TRAIT, FEATURE and CLASSIFIER (for the COMBINED score)
ResultType | Percentile, prevalence or likelihood
DataType | Detailed description of the category being evaluated
Value | The sample's value for this category if applicable
RefCancerType | Reference cancer type evaluated against
RefValue | Probability of the reference cancer type for this category of data

