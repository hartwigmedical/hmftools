# Relevant Oncogenic Summary of Eligibility
ROSE match all the genomic events and signatures that are determined by the Hartwig pipeline to the actionability treatment options in the Netherlands. 

## Contents
- [What is present in the actionability database?](#actionability-db)
- [How is actionability matched against genomic events and signatures?](#matching-of-actionability)
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

The field type_alteration could contain the following values:

Field  | Description
---|---
ONLY_HIGH | Only using actionable evidence for variants with high driver likelihood
ALWAYS | Always using actionable evidence
ALWAYS_NO_ACTIONABLE | Always add events to summary unless they aren't actionable events
OTHER | Other general messages to events

The field condition could contain the following values:

Field  | Description
---|---
AMPLIFICATION | The gene is actionable for amplification
EXTRACELLULAR_DOMAIN_MUTATION | ???
FUSION | The gene is actionable when the gene is present on 3' promiscuous or 5' promiscuous as fusion partner
INACTIVATION | The gene is actionable when the event inactivated the gene
INTERNAL_DELETION |
KINASE_DOMAIN_DUPLICATION |
LOSS |
POSITIVE |
RESISTANCE_MUTATION |
PURITY |
PURITY_UNRELIABLE |
FINDINGS |
GERMLINE | T
CUPPA |
CUPPA_INCONCLUSIVE |
NO_ONCOGENIC |
NO_ACTIONABLE |
NO_HRD_CAUSE |
NO_MSI_HRD_PROFILE |
NOT_BIALLELIC |

## Matching of actionability

Genomic events are categorized in six categories and actionability is matched for every category independently.

#### SNVs and (small) INDELs

For small variants (SNVs and INDELs) determined by [PURPLE](../purple/README.md) the type alteration will be annotated as follows:
- If the small variant is present in an ONCO gene the type alteration will be ACTIVATING_MUTATION
- If the small variant is present in an TSG gene the type alteration will be INACTIVATION

The actionability will be match, when meets with the following criteria: 
- The condition is ONLY HIGH for ACTIVATING_MUTATION or INACTIVATION the actionability will only be match for high driver variants
- The condition is ALWAYS_NO_ACTIONABLE for INACTIVATION the actionability will be match for all small variant unless the driver likelihood 

Do note that germline and somatic variants are treated equally. It is not considered relevant for clinical evidence whether the variant is
present in the germline already or has been acquired by the tumor somatically.

#### Copy numbers

Actionability on amplifications and deletions is considered applicable in case a gene has been classified as amplified or deleted by
[PURPLE](../purple/README.md). If we called a CNV the type alteration will be annotated as follows: 
- If it is an amplification the type alteration will be AMPLIFICATION
- If it is a deletion the type alteration will be LOSS

The actionability will be match, when meets with the following criteria:
- The condition is ALWAYS for AMPLIFICATION or LOSS the actionability will be always match 

#### Homozygous disruptions

When a gene has been homozygously disrupted according to [LINX](../linx/README.md), the type alteration will be INACTIVATION. 

The actionability will be match, when meets with the following criteria:
- The condition is ALWAYS for INACTIVATION for matching the actionability 

#### Fusions

For fusions that are deemed reportable according to [LINX](../linx/README.md) the following matching is performed:
- Fusions will be match to the actionability unless the driver likelihood 
- If fusion type is EXON_DEL_DUP this will match to INTERNAL_DELETION
- If fusion type is EXON_DEL_DUP and gene is EGFR with the specific range (25;26 exon up and 14;18 exon down) it will match to KINASE_DOMAIN_DUPLICATION
- If fusion type is PROMISCUOUS_3, PROMISCUOUS_5, KNOWN_PAIR, IG_KNOWN_PAIR or IG_PROMISCUOUS this will match to FUSION

#### Viral presence

For matching viral presence to actionability, the interpretation by [Virus Interpreter](../virus-interpreter/README.md) is used. 
The type alteration will be POSITIVE and the virus should meet with the following criteria: 
- The actionability will be always match because the condition is ALWAYS

#### Signatures
The signatures are categorized in four categories and actionability is matched for every category independently.

###### HRD
When a tumor has the signature homologous recombination repair (>= 0.5) according to [CHORD](https://github.com/UMCUGenetics/CHORD) 
the signature is match for actionability. 

###### MSI
When a tumor has the signature microsatellite instability (>= 4) according to [PURPLE](../purple/README.md) the signature is match 
for actionability.

###### TML
When a tumor has the signature tumor mutational load (>= 140) according to [PURPLE](../purple/README.md) the signature is match
for actionability.

###### TMB
When a tumor has the signature tumor mutational burden (>= 10) according to [PURPLE](../purple/README.md) the signature is match
for actionability.

#### Other
Next to 

Evidence on signatures is matched based on the comparator and cutoff defined by the evidence rule.
If the evidence rule provides no comparator and cutoff, the interpretation of the algorithm producing the signature is used to match.

## ROSE output

ROSE produces a tsv with the clinical eligibility of that sample

Field  | Description | Example
---|---|---
tumor sample | Sample ID                             | HMF123
conclusion | The clinical conclusion of the sample | BRCA2 inactivation, potential benefit from PARP inhibitors <enter> BRCA1 inactivation, potential benefit from PARP inhibitors <enter>

## Known issues 
- When we called two or more actionable variants, this shows 2 entries in the conclusion rather than merge this to one single sentence 
- We don't interpret the type alteration for NO_MSI_HRD_PROFILE, EXTRACELLULAR_DOMAIN_MUTATION and RESISTANCE_MUTATION

## Version History and Download Links
- [Upcoming]
    - Initial release. 