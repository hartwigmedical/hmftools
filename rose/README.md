# Relevant Oncogenic Summary of Eligibility
ROSE match all the genomic events and signatures that are determined by the Hartwig pipeline to the actionability treatment options in the Netherlands. 

## Contents
- [What is present in the actionability database?](#actionability-database)
- [How is actionability matched against genomic events and signatures?](#matching-of-actionability-based-on-genomic-events-and-signatures)
- [What output is produced by ROSE?](#rose-output)
- [Which known issues are present in ROSE?](#known-issues)
- [Version history and download links](#version-history-and-download-links)

## Actionability database
For determining which genomic event or signature is actionable in the NL an actionability database is used. This database has the following format:

Field  | Description | Example
---|---|---
match | The gene name or signature which could be actionable together with the type alteration | RAD51B
type_alteration | The alteration of the event |  INACTIVATION
condition | The condition where the event should meet with | ONLY_HIGH
conclusion | The clinical sentence which will be added into the conclusion | inactivation, potential benefit from PARP inhibitors (clinical trial)

The field condition could contain the following values:

Field  | Description
---|---
ONLY_HIGH | Only using actionable evidence for genomic events with high driver likelihood
ALWAYS | Always using actionable evidence unless driver likelihood 
ALWAYS_NO_ACTIONABLE | Always add genomic events to summary unless they aren't actionable events
OTHER | Other general messages to events

The field type alteration could contain the following values which are specific for genomic alterations/alterations:

Field  | Description
---|---
AMPLIFICATION | Actionable when the gene is amplified 
EXTRACELLULAR_DOMAIN_MUTATION | Actionable when a mutation is in the extracellular domain domain 
FUSION | Actionable when the fusion is KNOWN_PAIR, PROMISCUOUS_3, PROMISCUOUS_5, IG_KNOWN_PAIR or IG_PROMISCUOUS as annotated by PURPLE
INACTIVATION | Actionable when gene is inactivated 
INTERNAL_DELETION | Actionable when there is a EXON_DEL_DUP fusion of that gene 
KINASE_DOMAIN_DUPLICATION | Actionable when the genomic event determined a kinase domain duplication
LOSS | Actionable when the gene is deleted
POSITIVE | Actionable when there is a positive call for a specific signature e.g. HRD deficient 
RESISTANCE_MUTATION | Actionable when it is related to a resistence mutation

The field type alteration could contain also the following values which are related to disclaimers or general information which is 
important to know for interpretation of the clinical relevance:

Field  | Description
---|---
PURITY | Disclaimer sentence when there is a lower tumor purity low (below the 20%) 
PURITY_UNRELIABLE | Disclaimer sentence when the tumor purity could be relaible determined 
FINDINGS | A general sentence for what can be detected in the report
GERMLINE | A sentence if the small variant is detected in the germline of the patient but also called in the tumor 
CUPPA | The molecular tissue of origin classifier prediction with likelihood >= 80%
CUPPA_INCONCLUSIVE | The molecular tissue of origin classifier prediction which couldn't determined with enough likelihood 
NO_ONCOGENIC | The WGS analyse could not detected any oncogenic genomic event/signature
NO_ACTIONABLE | The WGS analyse could not detected any oncogenic genomic event/signature wich is actionable 
NO_HRD_CAUSE | The WGS analys detected an HRD mutational profile but no support is found for this signature 
NO_MSI_HRD_PROFILE | The WGS analys detected HRD/MSI support but the mutational profile isn't detected 
NOT_BIALLELIC | The WGS analyse called a TSG small variant but the variant itself isn't biallelic  

## Matching of actionability based on genomic events and signatures
Genomic events and signatures are categorized in six categories and actionability is matched for every category independently.

#### SNVs and (small) INDELs
Small variants (SNVs and INDELs) are determined by [PURPLE](../purple/README.md). For the small variants, the variants should be meets 
with the following conditions:
- For variants which are present in ONCO genes: 
  - The type alteration will be annotated with ACTIVATING_MUTATION
  - Actionability is only match when the variant is a high driver variant 

- For variants which are present in TSG genes:
  - The type alteration will be annotated with INACTIVATION 
  - Actionability is match when the variant is a high driver variant and condition is annotated with ONLY HIGH. Also, for some of called 
variants the precense is relevant for clinical interpretation but without actionability, the condition is annotated  ALWAYS_NO_ACTIONABLE
and will be added into the conclusion unless the driver likelihood. 
  - Checking the presence if the variant is bi-allelic. When an intact allele is still present, we notify this for interpretation

Do note that germline and somatic variants are treated equally. It is not considered relevant for clinical evidence whether the variant is
present in the germline already or has been acquired by the tumor somatically.

#### Copy numbers
Actionability on amplifications and deletions is considered applicable in case a gene has been classified as amplified or deleted by
[PURPLE](../purple/README.md). If an CNV is called, the actionability will be match with the following type alteration: 
- If it is an amplification the type alteration will be AMPLIFICATION
- If it is a deletion the type alteration will be LOSS

#### Homozygous disruptions
When a gene has been homozygously disrupted according to [LINX](../linx/README.md), the actionability will always match with the
type alteration INACTIVATION.

#### Fusions
For fusions that are deemed reportable according to [LINX](../linx/README.md), the following matching is performed:
- Fusions will be match to the actionability unless the driver likelihood 
- If fusion type is EXON_DEL_DUP this will match to the type alteration INTERNAL_DELETION actionability
- If fusion type is EXON_DEL_DUP and gene is EGFR with the specific range (25;26 exon up and 14;18 exon down) it will match to the 
actionability of KINASE_DOMAIN_DUPLICATION for EGFR
- If fusion type is PROMISCUOUS_3, PROMISCUOUS_5, KNOWN_PAIR, IG_KNOWN_PAIR or IG_PROMISCUOUS this will match to the type alteration
FUSION for actionability when the fusion gene is either on 3' promiscuous or 5' side 

#### Viral presence
For matching viral presence to actionability, the interpretation by [Virus Interpreter](../virus-interpreter/README.md) is used. 
The type alteration will be POSITIVE and the virus should be present in the actionability database. Also, the actionability is meets 
unless the driver likelihood

#### Signatures
The signatures are categorized in four categories and actionability is matched for every category independently.

###### HRD
When a tumor has the signature homologous recombination repair, which means a value >= 0.5 with the status omologous recombination 
deficient, according to [CHORD](https://github.com/UMCUGenetics/CHORD) the signature is match for actionability. Also, when there is no 
support for this signature this will be also mentioned. 

###### MSI
When a tumor has the signature microsatellite instability, which means a value >= 4 with the status MSI,  according to 
[PURPLE](../purple/README.md) the signature is match for actionability.

###### TML
When a tumor has the signature tumor mutational load, which means a value >= 140 with the TML high status, according to 
[PURPLE](../purple/README.md) the signature is match for actionability.

###### TMB
When a tumor has the signature tumor mutational burden, which means a value >= 10, according to [PURPLE](../purple/README.md) 
the signature is match for actionability.

#### Other matchings for the clinical conclusion

###### Molecular Tissue of Origin classifier
The molecular tissue of origin classifier predicts the primary tumor location of the patient based on the WGS date by 
[CUPPA](../cuppa/README.md). When this classify a tumor location with a likelihood above the 80% the primary tumor location is added. 
However, it is also possible that the likelihood is below the 80% and then 'inconclusive' is added to the clincial conclusion. 

###### Disclaimer of tumor purity
For every sample which we have analysed the tumor purity will be determined. For this, there are two flavours: 
- when the tumor purity couldn't be reliable determined, this will be mentioned
- when the tumor purity is below the 20% a disclaimer is added that the result should be interpret with caution

###### General information
In every clinical conclusion a sentence is mentioned that an overview of oncogenic DNA aberations can be found in the report. Next to this 
we could append some extra information: 
- When there are oncogenic DNA aberrations could be detected but those aberrations didn't actionable this will be mentioned in the summary
- When no oncogenic DNA aberrations could be detected this will be mentioned in the summary

## ROSE output

ROSE produces a tsv with the clinical eligibility of that sample

Field  | Description | Example
---|---|---
tumor sample | Sample ID | HMF123
conclusion | The clinical conclusion of the sample | BRCA2 inactivation, potential benefit from PARP inhibitors <enter> BRCA1 inactivation, potential benefit from PARP inhibitors <enter>

## Known issues 
- When we called two or more actionable variants, this shows 2 entries in the conclusion rather than merge this to one single sentence 
- We don't interpret the type alteration for NO_MSI_HRD_PROFILE, EXTRACELLULAR_DOMAIN_MUTATION and RESISTANCE_MUTATION

## Version History and Download Links
- [1.0](https://github.com/hartwigmedical/hmftools/releases/tag/rose-v1.0)
    - Initial release. 