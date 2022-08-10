# Search External Resources for Variant Evidence

SERVE harmonizes various sources of evidence into a single unified model that can be readily used in genomic analyses:  
 - A model is generated which allows mapping of genomic events to clinical evidence. 
 - An overview of mutations that are implied to be potential cancer drivers is generated.  

## Contents

* [Which knowledgebases are supported?](#supported-input-knowledgebases)
* [What is generated as final SERVE output?](#outputs)
* [How are genomic events extracted from the source knowledgebases?](#extraction-of-genomic-events-and-tumor-characteristics-from-knowledgebases)
* [What is done in terms of curation and harmonization?](#curation-and-harmonization-of-individual-knowledgebases)
* [What is done in terms of curation and harmonization of relevant treatment approaches?](#Relevant-treatment-approaches-of-the-evidence)
* [How does SERVE deal with multiple reference genome versions?](#handling-of-multiple-reference-genome-versions)
* [How does everything come together?](#overview-of-the-serve-algorithm)
* [Version history and download links](#version-history-and-download-links)

## Supported input knowledgebases

SERVE supports the ingestion of the following knowledgebases:
 - [CGI](https://www.cancergenomeinterpreter.org) - general purpose knowledgebase that is supported through [VICC](http://cancervariants.org)
 - [CIViC](https://civicdb.org) - general purpose knowledgebase that is supported through [VICC](http://cancervariants.org)
 - [CKB CORE](https://ckb.jax.org) - part of CKB's knowledgebase that is supported through [VICC](http://cancervariants.org)
 - [CKB FLEX](https://ckbhome.jax.org) - The complete CKB clinical database.
 - [OncoKB](https://www.oncokb.org) - general purpose knowledgebase that is supported through [VICC](http://cancervariants.org)
 - [DoCM](http://www.docm.info) - database containing pathogenic mutations in cancer
 - [iClusion](https://iclusion.org) - a database with all actively recruiting clinical trials in the Netherlands 
 - [ACTIN](https://github.com/hartwigmedical/actin/blob/master/treatment/README.md) - a database with all actively recruiting 
 clinical trials in the ACTIN study along with molecular inclusion criteria for these trials.
 - HMF Cohort - a database of recurrent somatic mutations in cancer-related genes from the Hartwig database.
 - HMF Curated - a database of known driver mutations curated by the Hartwig team.
 
Support for the following knowledgebases is under development:
 - [CBG Compassionate Use](https://www.cbg-meb.nl/onderwerpen/hv-compassionate-use-programma/overzicht-goedgekeurde-cup) - 
    a database of approved compassionate use programs in the Netherlands 
 
A number of other Hartwig modules support the ingestion (and analysis) of these knowledgebases:
 - [VICC Importer](../vicc-importer/README.md): A module supporting the ingestion of any knowledgebase ingested into VICC.
 - [iClusion Importer](../iclusion-importer/README.md): A client implementation of the iClusion API which transforms iClusion API output to 
 data that can be ingested into SERVE.
 - [CKB Importer](../ckb-importer/README.md): A module supporting the ingestion and analysis of CKB FLEX.  
  
Do note that SERVE does not provide the actual data that is input to the algorithm, but only provides support for its ingestion.
While SERVE itself is open-source, the sources that can be ingested have their own licensing and up to the users to make sure 
they are compliant with the usage of the data itself. 

## Outputs

SERVE generates clinical evidence in the following datamodel:
 - Treatment (name of trial or drug(s))
 - Relevant treatment approaches of the knowledgebase (the drug classes of the treatment related to the event)
 - Curated treatment approaches of the knowledgebase (the curated drug classes of the treatment related to the event)
 - Cancer type (annotated with DOID) for which the treatment is considered on-label.
 - Blacklist cancer types (annotated with DOID) that should be children of the main cancer type and are used for blacklisting 
 specific types of the main cancer type.
 - Tier / Evidence level of the treatment
 - Direction (Responsive for the treatment, resistant to the treatment or whether mutation implies no benefit for the treatment)
 - A set of URLs pointing towards the source website which provide extra information about the treatment.
 - A set of URLs with extra information about the evidence (e.g. publications backing up the evidence)
 
The following genomic events and tumor characteristics can be mapped to clinical evidence:
 - Genome-wide tumor characteristics such as signatures, MSI status, TML status or viral presence
 - Presence of specific HLA alleles
 - Multi-gene events such as gene fusions
 - Single gene events such as amplification or general (in)activation of a gene
 - Types of mutations in ranges overlapping with specific genes such as:
    - inframe insertions in EGFR exon 20
    - splice site mutations in MET exon 14
    - any type of missense mutation in BRAF codon 600
 - Specific missense mutations such as BRAF V600E
 
In addition to generating a mapping from various genomic events to clinical evidence, SERVE generates the following outputs describing 
genomic events implied to be able to driver cancer:
 - Specific known pathogenic fusion pairs
 - Known pathogenic amplifications and deletions
 - Known pathogenic exons (exons for which specific mutations are implied to be pathogenic)  
 - Known pathogenic codons (codons for which generic mutations are implied to be pathogenic)
 - Known pathogenic hotspots (specific mutations on specific loci)

## Extraction of genomic events and tumor characteristics from knowledgebases
 
### Gene checking

Evidence that is defined on a gene-level is checked to make sure that the gene exists in Hartwig's definition of the exome. 
If a gene does not exist in Hartwig's exome definition the evidence is ignored. For more information about Hartwig's definition 
of the exome, see [HMF Gene Utils](../gene-utils/README.md).

For fusions, genes are permitted that can exist in the context of a fusion pair (eg @IG genes). 

### Driver inconsistencies
In the supported knowledgebases, there can be events defined on genes that are are not part of the Hartwig's driver gene panel.
Also, the event could be inconsistent with respect to driver gene panel (e.g. "Inactivation" evidence for a gene that is configured to be 
an oncogene). There are various ways to deal with such inconsistencies on a per-knowledgebase level:

Filter  | Description
---|---
FILTER  | We filter every entry when the gene/event isn't present or there is an inconsistency with the Hartwig's driver gene panel
IGNORE  | Every gene/event is used regardless of mismatch/inconsistencies 
WARN_ONLY  | Every gene/event is used regardless of mismatch/inconsistencies, however a warning message is shown for the inconsistencies

### Protein resolving for SNVs and (small) INDELs
 
Evidence on SNVs and small INDELs generally come in their protein annotated form (e.g. BRAF V600E). 
SERVE uses [transvar](https://github.com/zwdzwd/transvar) to resolve these annotations into genomic coordinates (referred to as hotspots) 
for the reference genome version that is used by the input knowledgebase.
 
The first step is to choose what ensembl transcript to use for converting protein annotation back to genomic coordinates:
 1. If the knowledgebase configured a transcript for a mutation, that transcript is used exclusively.
 1. If no transcript is configured, SERVE uses the typical transcript used by Hartwig which is generally the canonical 
 transcript defined by ensembl.
 1. If a protein annotation does not exist on the canonical transcript and has no transcript configured in the knowledgebase, 
 a consistently specific transcript is picked for protein annotation in case multiple transcripts imply the same hotspot.
 
If a protein annotated form does not exist on any transcript for a specific gene, the evidence is ignored 
(see also [curation](#curation-and-harmonization-of-individual-knowledgebases)). 
 
Assuming a suitable transcript has been found, N hotspots are derived for each protein annotation as follows:
 - In case the mutation is caused by SNV or MNV every possible trinucleotide combination that codes for the new amino acid is generated.
 - In case the mutation is caused by a duplication (DUP) or an inframe deletion (DEL) 1 hotspot is generated which assumes 
 the exact reference sequence has been duplicated or deleted. 
   - In case a DEL can be left-aligned a hotspot is generated for every position between the left-aligned position and the actual position. 
 - In case the mutation is caused by an inframe insertion (INS) there are two flavors based on the length of the insertion:
   1. In case 1 amino acid is inserted, hotspots are generated for every trinucleotide coding for that amino acid.
   1. In case multiple amino acids are inserted, one of the potentially many hotspots is generated. This is just for practical reasons
    to put a limit on the (exponential) number of variants that can code for a multi-amino-acid insert. 
 - In case of a complex deletion/insertion (DELINS) the rules for hotspot generation for deletions and insertions are extrapolated. 
 Hence, the reference sequence is assumed to be deleted, and one new nucleotide sequence is inserted unless the insertion is 1 amino acid 
 in which case hotspots are generated for all trinucleotides coding for the inserted amino acid. 
 Complexity of the resulting variant is reduced by removing any bases that are shared between ref and alt at start or end of the variant.
 - In case of a frameshift the following hotspots are generated:
   - Any of the 12 possible single base inserts inside the affected codon that does not lead to synonymous impact in the affected codon
   - Any of the 3 possible single base deletes inside the affected codon that does not lead to synonymous impact in the affected codon
   - Any of the 2 possible double base deletes inside the affected codon that does not lead to synonymous impact in the affected codon
   
Additionally, hotspot generation is ignored for any INDEL that spans multiple exons. Examples are:
 - A DUP which duplicates a codon that is encoded by parts of two separate exons.
 - A frameshift which shifts into the intronic space of the gene.
  
Finally, Any INDEL longer than 50 bases is ignored since this is considered to be a structural variant rather than a small INDEL.
 
### Coordinate and mutation filter resolving for codons and exons
 
For evidence defined on a codon or exon level, no protein annotation resolving is done. 
Instead, genomic coordinates are resolved using the following rules. 

First off, evidence on codons and exons are assumed to be defined with respect to the Hartwig canonical transcript.  
  - If evidence for a specific codon or exon range is defined for a different transcript, this evidence is ignored. 
  Since all variants in Hartwig are annotated in terms of their impact on the Hartwig canonical transcript, resolving 
  this evidence could potentially lead to wrong matching.  
  - If no transcript is configured in the knowledgebase, it is assumed the canonical transcript is implied. 
 
For ranges that represent exons, the range is extended by 10 bases on both sides of the exon to be able to capture splice variants affecting 
the exon. 

In addition to resolving coordinates, every codon and exon range is annotated with a filter indicating which type(s) of mutations are valid 
for this range. SERVE tries to determine this based on the information specified in the knowledgebase, but if that information is not
sufficient, the Hartwig driver catalog is used to determine the filter.
  
Filter  | Description
---|---
NONSENSE_OR_FRAMESHIFT  | Only frameshifts or nonsense mutations are valid for this range
SPLICE | Only splice mutations are valid for this range
INFRAME | Any inframe INDEL (insert or delete) is valid for this range
INFRAME_DELETION | Only inframe deletions are valid for this range
INFRAME_INSERTION | Only inframe insertions are valid for this range
MISSENSE | Only missense mutations are valid for this range
ANY | Any mutation is considered valid for this range.

### Gene level event determination

For evidence that is applicable when a gene-wide level event has happened, the type of event required to match evidence to a 
mutation is derived from the knowledgebase event and no further interpretation is done. In case a knowledgebase provides insufficient 
details to make a decision (eg. mutant), we annotate to ANY_MUTATION. 

Gene level event  | Description
---|---
AMPLIFICATION  | Evidence is applicable when the gene has been amplified.
OVEREXPRESSION | Evidence is applicable when the gene has been amplified.
DELETION | Evidence is applicable when the gene has been completely deleted from the genome.
UNDEREXPRESSION | Evidence is applicable when the gene has been completely deleted from the genome.
ACTIVATION | Evidence is applicable when a gene has been activated. Downstream algorithms are expected to interpret this.
INACTIVATION | Evidence is applicable when a gene has been inactivated. Downstream algorithms are expected to interpret this.
ANY_MUTATION | SERVE does not restrict this evidence based on the type of mutation and considers every type of mutation applicable for this evidence.
FUSION | Evidence is applicable in case the gene has fused with another gene (either 3' or 5').
WILD_TYPE | Evidence is applicable in case no genomic alteration is detected.

### Exonic ranges specific for fusion pairs

For evidence on fusion pairs, SERVE can add restrictions on which exons are allowed to be involved in the fusion.
This is to support evidence on fusions like EGFRvII. 

Evidence on fusion pairs where these restrictions are missing can be assumed to be valid for any fusion between the two genes specified. 

### Genome wide tumor characteristics

For evidence that is applicable when a genome wide event has happened, the type of event required to match evidence to the event 
is derived from the knowledgebase event. When the knowledgebase event has a cutoff defined for this evidence this information will be also extracted. 
When no cut-off values is present but is expected for the characteristics, Hartwig's default cutoff values are used. 

Genome wide event  | Description
---|---
MICROSATELLITE_UNSTABLE  | Evidence is applicable when the genome has a MSI status
MICROSATELLITE_STABLE  | Evidence is applicable when the genome does not have a MSI status
HIGH_TUMOR_MUTATIONAL_LOAD | Evidence is applicable when the genome has a high tumor mutational load status
LOW_TUMOR_MUTATIONAL_LOAD | Evidence is applicable when the genome does not have a high tumor mutational load status
HIGH_TUMOR_MUTATIONAL_BURDEN | Evidence is applicable when the genome has a high tumor mutational burden status
LOW_TUMOR_MUTATIONAL_BURDEN | Evidence is applicable when the genome does not have a high tumor mutational burden status
HOMOLOGOUS_RECOMBINATION_DEFICIENT | Evidence is applicable when the genome has a HRD status
HPV_POSITIVE | Evidence is applicable when viral presence of some form of HPV has been found
EBV_POSITIVE | Evidence is applicable when viral presence of some form of EBV has been found

### Presence of HLA alleles 

For evidence on HLA Class Type I, no further interpretation is done. Only the exact HLA type should be known in 2 digits.

## Curation and harmonization of individual knowledgebases

Per knowledgebase curation and filtering is applied to harmonize knowledge from different sources and to correct/remove mistakes or 
evidence that is inconsistent with the Hartwig driver model.

### VICC Curation

For VICC the following curation and filtering is applied prior to presenting the data to SERVE:
 1. General filtering of mutations that are undetectable when analyzing DNA or RNA. Examples are phosphorylation and methylation.
 1. Determining whether the evidence is supportive of the specified direction. Eg if evidence "does not support" sensitivity
 we do not generate actionable results from this evidence.
 1. Filtering of specific mutations:
    - Mutations that remove the stop codon. These are simply not interpreted yet by the SERVE main algorithm.
    - Synonymous mutations in coding regions are assumed to be benign by SERVE and ignored.
    - Fusions that are not considered pathogenic by Hartwig are removed for lack of evidence of pathogenicity (regardless of their level of evidence).
    - Events that contradict Hartwig driver catalog. One example is "CCND3 loss" which is assumed to be benign. 
 1. Curation of specific mutations:
    - SNVs/INDELs that are not aligned correctly according to HGVS standards are corrected to be HGVS-compliant.
    - SNVs/INDELs that have correct notation but simply don't exist on the transcript specified by VICC are removed.
    - Fusion pairs for which the genes are in the wrong order are flipped around. 
    - Genes which are synonyms of genes used in the Hartwig exome definition are renamed.
 1. Correction of cancer types and DOID annotation:
    - Evidence for which DOID is missing have a DOID manually assigned.
    - Evidence on multiple cancer types generally get a wrong DOID assigned by VICC and are rectified.
 1. Correction of drugs for which A or B level evidence exists:
    - A whole range of drugs have wrong or inconsistent names in VICC and are rectified by SERVE.
    - VICC does not explicitly model the difference between "multiple different drugs" and a "combination treatment of multiple drugs". 
    This gets rectified by SERVE.

### DoCM Curation   

DoCM is used exclusively for known hotspot generation. The filtering is therefore tailored for hotspots:
 - Entries implying general codon mutations are removed.
 - Unusual notations for inframe deletions and insertions are removed.
 - Mutations that don't exist on the transcript specified by DoCM are removed.
 
Also, genes that do not follow HGNC model are renamed to their HGNC name.

### ACTIN Curation

ACTIN ingest the molecular inclusion criteria which are extracted from the ACTIN treatment database (see also [actin](https://github.com/hartwigmedical/actin/blob/master/serve-bridge/README.md)).
The inclusion criteria in the trials of the ACTIN database are defined in terms of specific rules.

Rule | When does a patient pass evaluation?
---|---
ACTIVATION_OR_AMPLIFICATION_OF_GENE_X | Activating mutation or amplification is found in gene X
ACTIVATING_MUTATION_IN_GENE_X | Activating mutation is found in gene X
FUSION_IN_GENE_X | Driver fusion with fusion partner gene X is found
SPECIFIC_FUSION_OF_X_TO_Y | Driver fusion with 2 specified fusion partner genes is found
INACTIVATION_OF_GENE_X | Inactivating mutation or deletion/disruption is found in gene X
MUTATION_IN_GENE_X_OF_TYPE_Y | Specific mutation Y is found in gene X
AMPLIFICATION_OF_GENE_X | Amplification is found in gene X
DELETION_OF_GENE_X | Deletion is found in gene X
WILDTYPE_OF_GENE_X | No driver mutation is found in gene X
MSI_SIGNATURE | MS Status = MSI
HRD_SIGNATURE | HR Status = HRD
TMB_OF_AT_LEAST_X | Tumor Mutational Burden (TMB) should be => X
TML_OF_AT_LEAST_X | Tumor Mutational Load (TML) should be => X
TML_OF_AT_MOST_X | TML should be <= X
HAS_HLA_A_TYPE_X | Patient should have at least one HLA allele of type X   

SERVE configures every trial to B-level evidence with `RESPONSIVE` direction in case the rule is involved in inclusion, and `NO_BENEFIT` in
case the rule is used for exclusion. The filtering is predominantly configurable rather than fixed in SERVE. 
The following filters can be configured in ACTIN:

Filter  | Description
---|---
FILTER_RULE_ON_GENE | Can be used to remove evidence of a specific rule of a particular gene 
FILTER_MUTATION_ON_GENE | Can be used to remove evidence of a gene with a specific mutation
FILTER_EVERYTHING_FOR_GENE | Can be used to remove all evidence of a gene (eg. Mutations that are inconsistent with the Hartwig driver catalog)
FILTER_EVERYTHING_FOR_RULE | Can be used to remove all evidence of a specific rule 

### iClusion Curation

iClusion contributes to actionability only. SERVE configures every trial to B-level evidence with `RESPONSIVE` direction. 
SERVE only considers trials with one or more molecular inclusion criterium. The filtering is predominantly configurable rather than fixed 
in SERVE. The fixed curation of iClusion that is done in SERVE is mapping gene names and signatures. 

The following filters can be configured in iClusion:

Filter  | Description
---|---
FILTER_EVENT_WITH_KEYWORD  | Can be used to remove evidence of a type that is not observable in DNA (eg "hypermethylation")
FILTER_VARIANT_ON_GENE  | Can be used to remove evidence of a gene (eg. Mutations that are inconsistent with the Hartwig driver catalog)
 
Finally, cancer types for which no DOIDs have been specified get a DOID assigned by SERVE.

### CKB FLEX Curation

For CKB FLEX curation and filtering is predominantly configurable rather than fixed in SERVE. The only fixed curation done in SERVE
is mapping evidence for tumor characteristics (such as MSI or High TMB) to actual characteristics since CKB FLEX models this as "genes".

The following filters can be configured for CKB FLEX, along with an example of how this is used by Hartwig:
Filter  | Description
---|---
ALLOW_GENE_IN_FUSIONS_EXCLUSIVELY  | CKB FLEX uses a hierarchy of events in such a way that every "fusion" is a child of "mutant". For certain genes (eg @IG) we want to ignore the abstract level and only include the fusion evidence since we only handle @IG on a fusion level in the Hartwig pipeline. 
FILTER_EVENT_WITH_KEYWORD  | Can be used to remove evidence of a type that is not observable in DNA (eg "hypermethylation")
FILTER_EXACT_VARIANT_FULLNAME  | Any specific variant can be removed through this filter. This is primarily used to remove variants that have a coding impact on their configured refseq transcript in CKB but are non-coding or don't exist on Hartwig's ensembl transcript.
FILTER_ALL_EVIDENCE_ON_GENE  | Is primarily used to remove evidence on genes which are simply not modeled correctly in Hartwig's gene model and hence can't be mapped properly
FILTER_EVIDENCE_FOR_EXONS_ON_GENE  | Some genes may have evidence on specific exons which don't exist on the ensembl transcript used by Hartwig
FILTER_SECONDARY_GENE_WHEN_FUSION_LEG  | Usage of this filter is similar to the use case for removing all evidence on genes. 

## Relevant treatment approaches of the evidence
External knowledgebases can be annotate evidence (treatment/event) with the relevant treatment approach. For making this usuable downstream, this 
is curated for harmonize the knowledge. The following filters can be configured: 

Filter  | Description
---|---
TREATMENT_APPROACH_CURATION  | The treatment approaches which should curated for downstream usages
DIRECTION_TREATMENT_APPROACH_CURATION_IGNORE | The treatment approach wouldn't be used because we didn't use the treatment approach for this direction 
EVENT_TREATMENT_APPROACH_CURATION_IGNORE | The treatment approach wouldn't be used because event is ignored for further interpretation 

## Handling of multiple reference genome versions
 
External knowledgebases generally define their knowledge for one specific reference genome version (v37 or v38). SERVE merges knowledgebases 
defined in either v37 or v38 reference genome versions. In addition SERVE generates its output for both reference genome v37 and v38.

### Ref-genome dependent knowledge extraction

Knowledge is extracted from any knowledgebase with respect to the ref genome version of that knowledgebase. The following resources
are used in a ref-dependent manner during knowledge extraction:
 - The reference genome fasta file 
 - The definition of Hartwig driver genes
 - The definition of Hartwig known fusions
 - The Hartwig ensembl data cache 

### Ref-genome dependent output

Once per-knowledgebase extraction is done, the extraction results are merged into a v37 and v38 version.
Any extraction for a v37 knowledgebase is taken over unchanged in the v37 output, and the same holds for any v38 knowledgebase into the v38 
output. For the remaining cases, the following conversion algo is executed.

#### Genomic position lift-over

Hotspots and ranges are lifted over using [HTSJDK's implementation](https://gatk.broadinstitute.org/hc/en-us/articles/360037060932-LiftoverVcf-Picard-) of UCSC LiftOver.
In case lift-over could not be performed a warning is raised unless the position is known to not exist in the target ref genome. 
In addition, any lift-over that lifts a position towards a different chromosome is raised as a warning but is accepted nonetheless.

There are a few additional checks for specific types of knowledge:
 - A hotspot lift-over is  accepted only in case the reference at the lifted position has remained unchanged.
 - A codon lift-over is accepted only in case the lifted range has a length of 3 bases.
 - Transcripts and codon/exon indices are removed from known codons and exons since they can't be trusted anymore after lift-over
 - In case the genomic region of a gene has been flipped between v37 and v38 we exclude the gene from liftover. 
 
## Overview of the SERVE algorithm

Every knowledgebase can be enabled or disabled through configuration.
SERVE starts with reading the various knowledgebases which are enabled. Knowledge is extracted after applying filtering and curation.
A knowledgebase can contribute to known and/or actionable events. Current configuration as follows:
      
Knowledgebase  | Ref genome version | Contributes to known events? | Contributes to actionable events?
---|---|---|---
ACTIN | v37 | No | Yes
CKB FLEX | v38 | Yes | Yes
DoCM | v37 | Yes | No 
Hartwig Cohort | v37 | Yes | No
Hartwig Curated | v37 | Yes | No
iClusion | v37 | No | Yes
VICC | v37 | Yes | Yes

Knowledge extraction is performed on a per-knowledgebase level after which all events are consolidated as follows:
 - All known events are aggregated on a per-event level where every event has a set of knowledgebases in which the event has been defined as pathogenic.
 - All actionable events are concatenated. Every actionable event that is present in multiple knowledgebases will be present multiple times in 
 the actionable output. 
  
 Within the Hartwig pipeline, SERVE output is used in the following manner:
  - The known output is used in various algorithms for various purposes. For example, the known hotspots produced by SERVE are used by 
  [SAGE](../sage/README.md) as the definition of the highest tier of calling and HOTSPOT annotation.
  - The actionable output is the database that [PROTECT](../protect/README.md) bases its clinical evidence matching on.
  
## Version History and Download Links
- [1.12](https://github.com/hartwigmedical/hmftools/releases/tag/serve-v1.12)
  - Support for generating a mySQL database
  - Support for adding the relevant treatment approaches to the output (which contains only drug classes when the source is CKB)
  - Support the distinction of amplification/overexpression and deletion/underexpression
  - When the knowledgebase defined evidence related to mutant, this evidence isn't anymore interpret as inactivation/activated related to TSG/ONCO gene
- [1.11](https://github.com/hartwigmedical/hmftools/releases/tag/serve-v1.11)
  - Use correct blacklisted tumor locations for solid tumors 
  - Update bug in mutation type filter for exon insertions and deletions
  - Remove "negative" from list of CKB events which are interpreted as inactivation events.
  - Move TMB evidence from CKB to evidence for tumor mutational burden rather than mutational load.
  - Various updates are made to ingestion of ACTIN source
  - Various updates are made to ingestion of iClusion source
    - Hotspot ERBB2 P780-Y781insGSP will be interpreted correctly 
    - Interpret EXON XX LOSS as internal deletion 
    - Interpret EXON XX DELETION as exon mutation 
- [1.10](https://github.com/hartwigmedical/hmftools/releases/tag/serve-v1.10)
  - Solve issues of v1.9
- [1.9](https://github.com/hartwigmedical/hmftools/releases/tag/serve-v1.9)
  - Support the raw input string of the input knowledgebases to the actionable output files
  - Support wild-type events as gene level evidences 
  - Support HLA Class type I as new actionability options 
  - The filtering of iClusion events is moved to a input file instead inside SERVE  
  - Created a link of CKB of the evidence for CKB Boost (web based)
  - For actionable signatures evidences could be applicable with different cut-offs. Now supporting those cut-off values of the different signatures (eg. TML >= 140 )
  - Support the possibility to blacklist specific tumor locations for particular treatments
  - Add an option to filter evidences when there are driver inconsistencies 
  - Support for curation the coordinates of genes because with ensembl data cache BRAF has the wrong coordinates
  - We support the interpretation of the new cancer type(DOIDs) of the CKB knowledgebase (JAX:10000009 and JAX:10000008)
- [1.8](https://github.com/hartwigmedical/hmftools/releases/tag/serve-v1.8)
  - Add full support for HGNC gene model implied by [HMF Gene Utils](../gene-utils/README.md)
    - Removed gene name mapping from ref genome converter when mapping between 37 and 38 (since gene names are always equal)
    - Remove gene name mapping from CKB gene extractor since CKB follows HGNC as well.
    - The ensembl data cache is used for resolving genes and canonical transcripts
    - Transvar uses new temporarily "new to old" gene mapping
    - Support for gene mapping in DoCM which does not follow HGNC
  - Add range annotation to (actionable) range:
    - Rename exonIndex to exonRank in KnownExons
    - Rename codonIndex to codonRank in KnownCodons
    - Add transcript, rangeType and rank to ActionableRange
  - Support for HRD in CKB
- [1.7](https://github.com/hartwigmedical/hmftools/releases/tag/serve-v1.7)
  - Extend config file from source iClusion
  - Extend config file from source CKB
  - "Advanced solid tumor" in iClusion is mapped to DOID 162 rather than 0050686 to avoid missing it for tumors with unknown tumor type
- [1.6](https://github.com/hartwigmedical/hmftools/releases/tag/serve-v1.6)
  - Add filter in VICC extraction to ignore evidence that does not support the direction when generating actionability.
- [1.5](https://github.com/hartwigmedical/hmftools/releases/tag/serve-v1.5)
  - "Advanced solid tumor" in CKB is mapped to DOID 162 rather than 0050686 to avoid missing it for tumors with unknown tumor type
- [1.4](https://github.com/hartwigmedical/hmftools/releases/tag/serve-v1.4)
  - Various additional checks to ref genome lift-over (such as filtering of events on genes for which strand has flipped).
  - CKB FLEX filtering framework has been added.
  - Solve bug when generating hotspots from MNVs that cross exon boundaries
- [1.3](https://github.com/hartwigmedical/hmftools/releases/tag/serve-v1.3)
  - Support for merging sources that differ in ref genome version (v37 vs v38).
  - Support for generating output for both ref genome version v37 and v38.
  - Driver catalog warnings are disabled for VICC.
  - KnownExons and KnownCodons are sorted more explicitly to make sure files don't change upon identical input.
  - An index file is generated for KnownHotspots VCF (for both v37 and v38). 
- [1.2](https://github.com/hartwigmedical/hmftools/releases/tag/serve-v1.2) 
  - Consistently pick a specific transcript for hotspot annotation in case multiple transcripts imply the same hotspot.
  - Extend splice sites from 5 bases to 10 bases beyond exon boundaries.
  - Add support for evidence for actionable viral presence (starting with EBV and HPV presence).
  - Add support for evidence of absence of high TMB and MSI.
  - Renamed actionable signatures to actionable characteristics.
- [1.1](https://github.com/hartwigmedical/hmftools/releases/tag/serve-v1.1) 
  - Ability to switch every resource on/off independently.
  - More predictable sorting of knowledge to ensure identical output on identical input.
  - Ability to curate VICC evidence levels in case they are suspicious.
- [1.0](https://github.com/hartwigmedical/hmftools/releases/tag/serve-v1.0)
  - Initial release.