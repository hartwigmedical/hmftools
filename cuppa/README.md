# Cuppa

## Overview
CUPPA is a tool that weighs multiple features observed from WGS and/or WTS data to predict the tissue of origin of a tumor sample. It is intended to 1) provide molecular tumor type prediction to verify histopathological classification, 2) provide support for specific tumor type classification in case of inconclusive histopathological outcome (differential diagnosis) and 3) prediction of primary tumor location for Cancer of Unknown Primary (CUP). 

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

Otherwise they can be obtained from the HMF Resources page: [HMFTools-Resources](https://resources.hartwigmedicalfoundation.nl/).

Cuppa will attempt to load these files from a reference directory supplied in its command-line arguments. The files are:

Data Type | Source | Filename | Description 
---|---|---|---
Sample Set | ALL | cup_ref_sample_data | List of the reference samples - SampleId and CancerType
Features | DNA | cup_ref_feature_prev.csv | Prevalence of fusions, viral insertions, drivers and known INDELs per cancer type
Features | DNA | cup_ref_driver_avg.csv | Average driver counts per cancer types 
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
  -cohort_sample_traits_file cup_ref_cohort_traits_data.csv \
  -cohort_sig_contribs_file cup_ref_cohort_signature_data.csv \
  -cohort_sv_data_file cup_ref_cohort_sv_data.csv \
  -cohort_features_file cup_ref_cohort_feature_data.csv \ 
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

# Algorithm

## Determination of cancer type cohorts

Cohorts for training the algorithm were constructed from the HMF database by selecting the highest purity sample from each unique patient from our database with qcStatus = ‘PASS’. 37 tumor categories were defined based on the clinical annotations in the HMF database of primaryTumorLocation, primaryTumorSubLocation, primaryTumorType and primaryTumorSubType as follows:

CUPPA Category | primaryTumorLocation:subLocation (primaryTumorType:primaryTumorSubType)
---|---
Acute myeloid leukemia | 
Anogenital | Penis, Vulva, Vagina, Anus ({excl. melanoma}), Uterus:Cervix
Bile duct /Gallbladder | Bile duct; Hepatobiliary system; Gallbladder
Bone/Soft tissue: Other | Bone/Soft tissue ({other or unspecified})
Breast | Breast
Cartilaginous neoplasm | 
Chronic lymphocytic leukemia |
Colorectum/Appendix/Small intestine | Colorectum ({other or unspecified}); Appendix; Small intestine({other})
Esophagus/Stomach | Esophagus ({excl. Nueroendocrine tumor}); Stomach ({excl. Nueroendocrine tumor}); Gastroesophageal
GIST | Bone/Soft tissue (Gastrointestinal stromal tumor)
Glioma | Nervous system (Glioma)
Head and neck: other | Head and neck({other})
Kidney | Kidney
Leiomyosarcoma | Bone/Soft tissue (Leiomyosarcoma)
Liposarcoma | Bone/Soft tissue (Liposarcoma)
Liver | Liver ({excluding Nueroendocrine tumor})
Lung: NET | Lung(Neuroendocrine tumor)
Lung: Non-small Cell | Lung(Carcinoma:Non-small cell carcinoma); Lung(Carcinoma:Adenocarcinoma); Lung({other})
Lung: Small Cell | Lung(Carcinoma:Small cell carcinoma); Lung(Carcinoma:Small cell carcinoma combined type)
Lymphoid tissue | Lymphoid tissue
Medulloblastoma | Medulloblastoma
Mesothelium | Mesothelium
Osteosarcoma | Bone/Soft tissue (Osteosarcoma)
Other | Gastrointestinal tract, Eye, Bone marrow, Nervous system({other}), Adrenal Gland, Thymus, Testis, Esophagus (Nueroendocrine tumor),Stomach (Neuroendocrine tumor)
Ovary/Fallopian tube | Ovary; Fallopian tube
Pancreas | Pancreas ({other>)
Pancreas: NET | Pancreas (Neuroendocrine Tumor)
Pilocytic astrocytoma | Nervous system (Pilocytic astrocytoma)
Prostate | Prostate
Salivary gland/Adenoid cystic | Head and Neck:Salivary gland, Head and Neck:Parotid gland, Head and Neck:Sublingual gland, {any}(Carcinoma:Adenoid cystic carcinoma), Trachea
Melanoma | <Any, excluding Eye> (Melanoma)
Skin:Other | Skin ({other})
Small intestine/Colorectum: NET | Small intestine(Neuroendocrine tumor); Colorectum(Neuroendocrine tumor)
Thyroid gland | Thyroid gland
Urothelial tract | Urothelial tract
Uterus:Endometrium | Uterus:Endometrium
Myeloproliferative neoplasm |

Certain cancers such as Esophagus and Stomach were combined for the categorisation as we found empirically that the CUPPA classifiers had little ability to distinguish between them. For other cancers including Lung, Bone/Soft tissue, Skin, Uterus & Pancreatic cancers we have broken into subtypes where histological information allows. All cancers not in one of these 36 cohorts was deemed as “Other” and was excluded from the reference cohorts for analysis.   Samples with ‘unknown’ tumor type are also excluded. Finally, 45 samples were also explicitly excluded from the reference cohort where our analysis strongly suggested the clinical configured cancer type may be incorrect for these samples

## DNA Classifier logic

CUPPA includes 3 orthogonal DNA classifiers based on positional mutational distribution, SNV mutational profile and feature prevalence, and a 4th classifier which combines the 3 together to make an overall prediction. Each classifier assigns a likelihood to each cancer type with the sum of the likelihoods adding up to 1 across the cancer types. 

The algorithm for each of the classifiers is described below:

### GENOMIC_POSITION_SIMILARITY

This classifier solely relies on the mutational distribution of tumors of genomic position, which has been shown previously to have strong predictive power for tissue of origin (eg. https://www.nature.com/articles/s41467-019-13825-8). 

CUPPA calculates a consensus mutation distribution for each cohort by counting SNV TMB by bucketed genomic position across each cohort. High TMB samples are downsampled to 20k mutations in this consensus so that individual samples cannot dominate a cohort. CUPPA counts mutations using a window size of 500kb bases (chosen after testing various sizes from 100kb to 10Mb). 

The genomic position similarity likelihood for a given sample is determined by first calculating the cosine similarity (CSS) of a sample to each cohort consensus distribution and then weighing using the following algorithm:
```
Score(sample=s,cancerType=i) = 8^[100*(CSS(i,s)-BestCSS(s))] 
```
CUPPA sums the scores across each tumor type to estimate a likelihood for each cancer type:
```
Likelihood(tumorType=i) = Score(i) / SUM(all tumors) [Score]
```

### SNV_96_PAIRWISE_SIMILARITY

This classifier relies solely on relative SNV counts via the 96 trinucleotide buckets frequently used for cosmic signatures. The cosmic signatures are not used directly, but the classifier is designed to capture the obvious similarities that can also be observed via signatures capturing known cancer specific mutagenic effects such as UV & Smoking and also background signatures per cancer type. 

Unlike the genomic position similarity which determines a consensus view of mutational distribution, the SNV_96_PAIRWISE classifier does not create a consensus view per tumor type as tumor types may have a diverse range of mutational profiles. Instead the classifier calculates a pairwise cosine similarity between the sample in question and every other sample in the Hartwig cohort.

Once a pairwise CSS has been determined, a score is calculated for each pair using the following formula:
```
Score(i,j) = 50^[-100*(1-CSS)] ^[ maxCSS^8] * mutationCountWeightFactor * cohortSizeWeightFactor
```
Where:
* MaxCSS is the maximum pairwise CSS for any sample in the cohort. This factor reduces confidences in general for samples that have no close pairwise match.
* mutationCountWeightFactor penalises pairs with large differences in SNV TMB. This is implemented as:
```
mutationCountWeightFactor = min(SNV_TMB(i)/SNV_TMB(j),SNV_TMB(j)/SNV_TMB(i)) 
```
* cohortSizeWeightFactor penalises larger cohorts which will have more similar tumors just by chance (eg. Breast cohort =~ 750 samples vs Thyroid cohort =~ 20 samples), implemented as:
```
cohortSizeWeightFactor = sqrt(# of samples of tumor type) / SUM(i)[sqrt(# of samples of tumor type i)]
```
As for genomic position similarity, CUPPA sums the scores across each tumor type to estimate the likelihood:
```
Likelihood(tumorType=i) = SUM(tumorType=i)[ Score] / SUM(all tumors) [Score]
```

### FEATURE 

The FEATURE classifier uses observed prevalence of both cancer type specific drivers as well as certain passenger mutational features that may be significantly enriched or depleted in certain types to predict the cancer type of a sample. 

#### Driver Prevalence
Driver (or driver like) features used include all driver point mutation, high amplification, homozygous deletion and homozygous disruptions in the driver catalog as well as viral insertions & fusions. For fusions, known pathogenic fusion pairs, IG rearrangement pairs and exon deletions/duplications configured in the HMF fusion knowledge base are all considered as features as are fusions with highly promiscuous exons such as ALK exon 20-21. For Sarcomas specifically, we override the prevalence for a list of 56 pathognomic fusions which are highly diagnostic but may not be prevalent enough to be present in our database to the appropriate cancer type with the maximal allowed feature weight.

Indels in repeat contexts of 6 or less bases in 3 lineage defining genes: ALB (highly specific to Liver cancer) and SFTPB & SLC34A2 (highly specific to Lung cancer) are also treated as additional features (note though that they are ignored for MSI samples). A set of Lung cancer specific EGFR hotspots (including T790M, L858R and exon 19 and 20 inframe deletions) are also treated as a single feature.

Features are weighted by driver likelihood. For point mutations the driver likelihood (the dnds calculated probability between 0 and 1 that the mutation is a driver) is used to weight the mutations, whilst other mutations, virus insertions and fusions are assumed to have probability of 1. 

The prevalence of each feature in each cancer type is calculated
```
Prevalence = minPrevalence + sum (driverLikelihood) / COUNT(samples)
```
Where minPrevalence is a fixed notional background rate of observing a passenger set to 0.15 / count of cancer types for drivers or indels in lineage defining genes and 0.01 / count of cancer types for fusions and viral insertions which are rarely passengers. 

A combined driver score for each cancer type is calculated by taking the product of the observed prevalence of each of the drivers from the sample in the cancer type cohort, discounted by the driver likelihood in the cancer itself. ie:
```
DriverScore = weightFactor(cohort)* PRODUCT[Prevalence(d)^driverLikelihood(d,s)]
```
Where the weight factor = meanDriverLoad(pan-cancer) / meanDriverLikelihood(cohort) and is intended to reduce the tendency for cancer types with higher average rates of drivers such as Urinary Tract and Esophagus to have higher driver scores

#### Passenger Prevalence
In addition to drivers, mutational burdens of certain types of events can vary widely across different cancer types. For example LINE insertions are universally observed in Esophagus and certain other cancers but almost non-existent in other cancers. Depending on the feature it may be useful to test that the rate observed is either higher or lower than what is expected of the cancer type. 

Since different cancers may have different characteristic frequencies, this is modeled for this classifier as a prevalence with a dynamic cutoff based on the rate observed in the sample itself. Specifically if testing for an enriched rate, the cutoff is set to 25% below the observed rate limited to a maximum value of the highest observed 95th percentile rate of any cancer cohort. Conversely if testing for a depleted rate, the cutoff is set to 25% below the observed rate limited to a maximum value of the highest observed 95th percentile rate of any cancer cohort. 

The following features are tested for enrichment and/or depletion:

Feature | Enrichment | Depletion
---|---|---
SNV_TMB | TRUE | TRUE
MS_INDEL_TMB | TRUE | TRUE
LINE_COUNT | TRUE | TRUE
TELOMERIC_SGL_BE_COUNT | TRUE | NA
MAX_COMPLEX_SIZE | TRUE | NA
SIMPLE_DUP_32B_200B | TRUE | NA

As for drivers the prevalence in each cancer type is added to a minPrevalence set to 0.15 / count of cancer types. The passenger score is simply the product of all the passenger prevalence rates
```
PassengerScore = PRODUCT[max(Passenger Prevalence,minPrevalence)]
```

#### Combining scores to a likelihood
The passenger and driver scores are multiplied together to get a single score:
```
Score = PassengerScore * Driver Score
```
And finally CUPPA sums the scores across each tumor type to estimate the likelihood:
```
Likelihood(tumorType=i) = Score(i)^correlationDampenFactor / SUM(all tumors) [Score^correlationDampenFactor]
```
The correlationDampenFactor is introduced to reduce the confidence of the classifier and set at 0.8 to empirically match the observed accuracy. This is required as some of the driver or passenger features may be correlated with each other - for example same arm amplifications are highly correlated and TMB might be positively correlated with more drivers in general

### DNA_COMBINED CLASSIFIER

A combined score is calculated by multiplying the 3 likelihoods together with an absolute floor set at 1% per likelihood. The likelihood is then calculated as
```
Likelihood(tumorType=i) = PRODUCT(max(0.01,Classifier(i,j)))^correlationDampenFactor / SUM(all tumors)[PRODUCT(max(0.01,Classifier(j)))]^correlationDampenFactor
```
As for the feature classifier, a correlationDampenFactor is introduced to reduce the confidence of the classifier and reflect the fact that the individual classifiers are not completely independent. A value of 0.65 is chosen to empirically match the confidence to the observed accuracy.

For the DNA_COMBINED classifier, males are excluded from matching ‘Ovary’ and ‘Uterus’ cancer cohorts and females are excluded from matching the ‘Prostate' cohort. ‘Breast’ cancer scores for male cancer cohorts are penalised but not excluded.

## RNA Classifier logic

CUPPA has 2 WTS based RNA classifiers and a combined RNA classifier:

### EXPRESSION_PAIRWISE classifier
The pairwise classifier calculates a pairwise cosine similarity of log(adjTPM+1) per gene, between the sample in question and every other sample in the Hartwig cohort.

Once a pairwise CSS has been determined, a score is calculated for each pair using the following formula:
```
Score(i,j) = 50^[-100*(1-CSS)] * cohortSizeWeightFactor
```
Where cohortSizeWeightFactor penalises larger cohorts which will have more similar tumors just by chance. It its calculated as:
```
cohortSizeWeightFactor = sqrt(# of samples of tumor type) / SUM(i)[sqrt(# of samples of tumor type i)]
```
CUPPA then sums the scores across each tumor type to estimate the likelihood:
```
Likelihood(tumorType=i) = SUM(tumorType=i)[Score] / SUM(all tumors)[Score]
```

### Novel Splice Junction (ALT_SJ_COHORT) classifier

A novel splice junction is defined in this context as any splice junction that is not annotated in ensembl.  A set of recurring novel splice junctions sites were identified within each cancer cohort - ie. those with 3 or more fragments supporting a novel site in 2 or more samples. A reference file was then formed by calculating the average fragment count per cancer cohort at each of these novel sites.

The novel splice junction classifier tests a sample’s fragment counts against each cancer cohort’s average fragment count per novel splice junction site. This is done by calculating a cosine similarity of log(fragmentCount + 1).
```
Score(sample=s,cancerType=i) = 3.5^[100*(CSS(i,s)-BestCSS(s))] 
```
CUPPA sums the scores across each tumor type to estimate a likelihood for each cancer type:
```
Likelihood(tumorType=i) = SUM(tumorType=i)[ Score] / SUM(all tumors) [Score]
```

#### Impact of read length
The HMF RNA cohort contains a mix of samples sequenced with 151 and 76 read lengths, and each of these lengths exhibit differences in novel splice junction fragment support. The 151-read-length samples were sequenced with greater depth, and in addition to often having greater fragment count support per novel splice site, approximately 10% of novel splice sites were only present in 151 read-length samples. Those cancer cohorts with a predominance of samples with either read-length of 76 or 151 read bases tended to find a closer CSS match with other samples of the same read-length. 

To address this bias, the reference cancer cohort file was split into average fragment counts per cancer type and per read-length group. A sample was only tested against the cancer reference data for its matching read length sub-cohorts.

### RNA_COMBINED classifier

A combined RNA classifier is calculated using the same formula as the combined DNA based on the 2 expression classifiers. The correlationDampenFactor is set to 0.7 via empirical analysis for the RNA_COMBINED confidence calculation. Gender restrictions are applied in the same manner as for the DNA_COMBINED score.

## Overall Combined Classifier

The RNA and DNA classifiers can be further merged into a consensus classifier in the same manner as the RNA_COMBINED and DNA_COMBINED, by merging all 5 individual classifiers. The correlationDampenFactor is set to 0.4 via empirical analysis for the overall COMBINED confidence calculation. Gender restrictions are applied in the same manner as for the DNA_COMBINED score.

## Nearest Neighbor Analysis

In addition to the classifiers, the 20 nearest neighbour samples by pairwise cosine similarity are reported for 3 different features:
* Count of SNV TMB per 500k genomic position buckets 
* Count of SNV TMB by 96 mutational context bucket
* Log(TPM+1) RNA expression by gene

Note that all samples are used for this analysis including rare cancer types that are not one of the CUPPA categorisations used in the classifiers.

## Known potential biases or issues

Bias | % of samples | Classifier | Description
---|---|---|---
AID_APOBEC | >2% | SNV_96 / GENOMIC_ POSITION  | Signature shared across 5-6 cohorts, but strongest in Urothelial / Breast. The genomic position signature for AID_APOBEC seems to be very different.  Lung and Eso/Stomach samples in particular get low GEN_POSITION. Other cancer types such as Anogenital & Head & neck perform ok on GEN_POS, but poorly on other classifiers.  
Small cohort size | 2% | All | Rounding issues and noise dominate all classifiers where cohort size is small (<25 samples), prevents us from small cohorts such as Testis, and diminishes performance even > 25 samples. Also true for pairwise classifiers even though we adjust for it.
TMBPerMb < .7 | 0.5% | ALL DNA | Generally low confidence.   Often mismatch to Pancreas:NET, likely due to ‘Low TMBPerMB’ feature
High driver load | 1% | FEATURE | Samples with a high number of drivers tend to match Urothelial Tract cancers (these have the highest rate of drivers)
MSI | 0.3% | GENOMIC_ POSITION | Samples with MSI typically have very low GENOMIC_POSITION scores to the correct cancer type.  Similar to AID_APOBEC effect 
Pathognomonic events | 1% | FEATURE | Rare pathognomonic events may not be found previously in our cohort or may not be weighed highly enough due to ‘min_prevalence’.  For some drivers the mechanism may be diagnostic whereas we only calculate features at a gene level, eg: SPOP amp (breast) vs mutation (prostate), KIT amp (lung) vs mutation (sarcoma), FOXA1 amp (lung) vs  mutation (breast/prostate), KRAS amp (esophagus) vs mutation (CRC/pancreas), Hypermutations in BCL2 & other genes (Lymphoid)
Metastasis site | 0.2% | RNA |  Liver mets can be mistaken for liver primary particularly low purity samples
Copy number | ? | GENOMIC_ POSITION | Adjusting for copy number may improve weightings
Treatment signatures | 0.1% | SNV_96 | Samples with strong treatment signatures (eg SYD985) will match each other with high certainty
Lung: Small-cell vs non small cell | 0.5% | GENOMIC_ POSITION | ‘Lung: non-small cell’ can strongly match the genomic position profile of Lung: small cell with high confidence.   Possibly due to timing of transformation to small cell?
Bile Duct / Gallbladder | 1% | ALL | Can be mistaken for Liver or Pancreas with high confidence
Non-smoking Lung | 0.2% | GENOMIC_ POSITION | Performance is weaker, but can mostly be explained by AID_APOBEC / pathognomonic events
Esophagus / Stomach vs Colorectal | 0.5% | ALL RNA | Esophagus frequently presents as Colorectal on all RNA classifiers
Anogenital vs Head & Neck: Other | 0.4% | All | Can often be mistaken for each other.
Sarcoma | 1.5% | ALL | Frequent mismatches between Leiomyosarcoma, Liposarcoma, Osteosarcoma and ‘other’.  Multiple causes.  Larger cohorts would help make clearer cohorts and could allow distinct groups for Rhabdomyosarcoma and others.  Some samples are marked as ‘Sarcoma’ and matched to Leiomyosarcoma are reported as match=F, but may be TP. Spindle cell sarcoma appear to group better with Leioymyosarcoma but are marked as ‘other’
Liposarcoma | ? | GENOMIC_ POSITION | MDM2+CDK4 coamplified liposarcomas (well-differentiated/dedifferentiated liposarcoma) resolve better to the Liposarcoma cohort compared to liposarcomas with diagnostic fusions (e.g. myxoid liposarcoma)



## Version History and Download Links
- [1.0](https://github.com/hartwigmedical/hmftools/releases/tag/cuppa-v1.0)
- [1.1](https://github.com/hartwigmedical/hmftools/releases/tag/cuppa-v1.1)
- [1.2](https://github.com/hartwigmedical/hmftools/releases/tag/cuppa-v1.2)
- [1.3](https://github.com/hartwigmedical/hmftools/releases/tag/cuppa-v1.3)
- [1.4](https://github.com/hartwigmedical/hmftools/releases/tag/cuppa-v1.4)
- [1.5](https://github.com/hartwigmedical/hmftools/releases/tag/cuppa-v1.5)
- [1.6](https://github.com/hartwigmedical/hmftools/releases/tag/cuppa-v1.6)







