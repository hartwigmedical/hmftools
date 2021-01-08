# **S**earch **E**xternal **R**esources for **V**ariant **E**vidence (SERVE)

SERVE harmonizes various sources of evidence into a single unified model that can be readily used in genomic analyses.  
The model provides a mapping from genomic events to clinical evidence. 
In addition, SERVE generates output containing all genomic events that are implied to be able to driver cancer.  

### Supported input knowledgebases

SERVE supports the ingestion of the following knowledgebases:
 - [CGI](https://www.cancergenomeinterpreter.org) - general purpose knowledgebase that is supported through [VICC](http://cancervariants.org)
 - [CIViC](https://civicdb.org) - general purpose knowledgebase that is supported through [VICC](http://cancervariants.org)
 - [CKB CORE](https://ckb.jax.org) - part of CKB's knowledgebase that is supported through [VICC](http://cancervariants.org)
 - [OncoKB](https://www.oncokb.org) - general purpose knowledgebase that is supported through [VICC](http://cancervariants.org)
 - [DoCM](http://www.docm.info) - database containing pathogenic mutations in cancer
 - [iClusion](https://iclusion.org) - a database with all actively recruiting clinical trials in the Netherlands 
 along with cancer types and molecular inclusion criteria for these trials.
 - HMF Cohort - a database of recurrent somatic mutations in cancer-related genes from the Hartwig database.
 - HMF Curated - a database of known driver mutations curated by the Hartwig team.
 
The following knowledgebases are under development:
 - [CBG Compassionate Use](https://www.cbg-meb.nl/onderwerpen/hv-compassionate-use-programma/overzicht-goedgekeurde-cup) - 
    a database of approved compassionate use programs in the Netherlands
 - [CKB FLEX](https://ckb.jax.org) - The complete CKB clinical database. 
 
A number of hmf modules support the ingestion (and analysis) of these knowledgebases:
 - [VICC Importer](../vicc-importer/README.md): A module supporting the ingestion of any knowledgebase ingested into VICC.
 - [iClusion Importer](../iclusion-importer/README.md): A client implementation of the iClusion API which writes the iClusion to a flat 
 file that can be ingested into SERVE.
 - [CKB Importer](../ckb-importer/README.md): A module supporting the ingestion (and analysis) of CKB FLEX.  
  
Do note that SERVE does not provide the actual data that is input to the algorithm, but only provides support for its ingestion.
While SERVE itself is open-source, the sources that can be ingested have their own licensing and up to the users to make sure 
they are compliant with the usage of the data itself. 

### Outputs

Clinical evidence is generated in the following datamodel:
 - Treatment
 - Cancer type (including DOID) for which the treatment is on-label.
 - Tier / Evidence level of the treatment
 - Direction (Responsive for the treatment or resistant to the treatment)
 
The following genomic events can be mapped to clinical evidence:
 - Genome-wide events such as signatures or MSI status
 - Multi-gene events such as gene fusions
 - Single gene events such as amplification or general (in)activation of a gene
 - Types of mutations in genic ranges such as:
    - inframe insertions in EGFR exon 20
    - splice site mutations in MET exon 14
    - any type of missense mutation in BRAF codon 600
 - Specific missense mutations such as BRAF V600E
 
In addition to generating a mapping from various genomic events to clinical evidence, SERVE also generates the following outputs describing 
genomic events implied to be able to driver cancer:
 - Specific known fusion genes including individual genes that are believed to be able to drive cancer when fused with any other gene 
 (in either 5' or 3' position)
 - Known amplifications and deletions
 - Known hotspots (specific mutations on specific loci)
 - Known codons (codons for which generic mutations are implied to be pathogenic)
 - Known exons (exons for which specific mutations are applied to be pathogenic)

SERVE can be configured to generate its output either for reference genome version 37 or version 38.  

## Extraction of genomic events from knowledgebases
 
### Gene checking

Any type of evidence that is ingested into SERVE has to be defined for a gene which is known in Hartwig's definition of the complete exome.

For fusions, genes are added that are only defined in the context of a fusion pair (eg IG genes) 
 
### Protein resolving for SNVs and (small) INDELs
 
Evidence on SNVs and small INDELs generally come in their protein annotated form (eg BRAF V600E). 
SERVE uses [transvar](https://github.com/zwdzwd/transvar) to resolve these annotations into genomic coordinates (referred to as hotspots) 
for the reference genome version that is configured to be used.
 
The first step is to choose what transcript to use for converting protein annotation back to genomic coordinates:
 1. If the knowledgebase configured a transcript for a mutation we exclusively use that transcript.
 1. If no transcript is configured, SERVE uses the typical transcript used by Hartwig which is generally the canonical 
 transcript defined by ensembl.
 1. If a protein annotation does not exist on the canonical transcript and has no transcript configured in the knowledgebase, 
 any random transcript for which the protein annotation does exist is picked.
 
Assuming we found a suitable transcript, we then derive N hotspots for each protein annotation as follows:
 - In case the mutation is caused by SNV or MNV every possible trinucleotide combination that codes for the new amino acid is generated.
 - In case the mutation is caused by a duplication (DUP) or an inframe deletion (DEL) 1 hotspot is generated which assumes 
 the exact reference sequence has been duplicated or deleted. 
   - In case a DEL can be left-aligned a hotspot is generated for every position between the left-aligned position and the actual position. 
 - In case the mutation is caused by an inframe insertion (INS) there are two flavors based on the length of the insertion:
   1. In case 1 amino acid is inserted, hotspots are generated for every trinucleotide coding for that amino acid.
   1. In case multiple amino acids are inserted, one of the potentially many hotspots is generated.
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
 
For resolving coordinates for codons and exons the Hartwig canonical transcript of a gene is used exclusively. If evidence for a specific codon
or exon range is defined for a different transcript this evidence is ignored. If no transcript is configured in the knowledgebase, it is 
assumed the canonical transcript is implied. 
 
For ranges that represent exons, the range is extended by 5 bases on both sides of the exon to be able to capture splice variants affecting 
the exon. 

In addition to resolving coordinates, every codon and exon range is annotated with a filter indicating which types of mutations are valid 
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
ANY | Any of the above mutations are valid for this range (synonymous mutations and other types of mutations are still excluded).

### Gene level event determination

For evidence that is applicable when something has happened on a gene level, the type of event that is required to match evidence to a 
mutation is derived from the knowledgebase event. In case a knowledgebase provides insufficient details to make a decision, the Hartwig 
driver catalog is used to determine what event qualifies for the evidence. 

Gene level event  | Description
---|---
AMPLIFICATION  | Evidence is applicable when the gene has been amplified.
DELETION | Evidence is applicable when the gene has been completely deleted from the genome.
ACTIVATION | Evidence is applicable when gene has been activated. Downstream algorithms are expected to further define this.
INACTIVATION | Evidence is applicable when gene has been inactivated. Downstream algorithms are expected to further define this.
ANY_MUTATION | This type is a "catch-all" meaning the knowledgebase did not provide enough details, and the gene is not a driver gene in Hartwig's driver catalog.
FUSION | Evidence is applicable in case the gene has fused with another gene.

### Exonic ranges specific for fusion pairs

For evidence on fusion pairs, SERVE can add restrictions on which exons are allowed to be involved in the fusion.
This is to support evidence on fusions like EGFRvII. 

Evidence on fusion pairs where these restrictions are missing can be assumed to be valid for any fusion between the two genes specified. 

## Curation and harmonization of individual knowledgebases

Per knowledgebase curation and filtering has been implemented to harmonize knowledge from different sources and to correct/remove mistakes.

### VICC Curation

VICC contributes to both actionable and known events.

For VICC the following curation and filtering is applied prior to presenting the data to SERVE:
 1. General filtering of mutations that are undetectable when analyzing DNA or RNA. Examples are phosphorylation and methylation.
 1. Filtering of specific mutations:
   - Mutations that remove the stop codon. These are simply not interpreted yet by the SERVE main algorithm.
   - Synonymous mutations in coding regions are assumed to be benign by SERVE and ignored.
   - Fusions that are not considered pathogenic by Hartwig are removed for lack of evidence of pathogenicity (regardless of their level of evidence).
   - Events that contradict Hartwig driver catalog. An example is "CCND3 loss" which is assumed to be benign. 
 1. Curation of specific mutations:
   - SNVs/INDELs that are not aligned correctly according to HGVS standards are corrected to be HGVS-compliant.
   - SNVs/INDELs that have correct notation but simply don't exist on the transcript specified by VICC are removed.
   - Fusion pairs for which the genes are in the wrong order are flipped around. 
   - Genes which are synonyms of genes used in the Hartwig exome definition are renamed.
 1. Correction of cancer types and DOID annotation:
   - Evidence for which DOID is missing is added manually.
   - Evidence on multiple cancer types generally get a wrong DOID assigned by VICC and these are rectified.
 1. Correction of drugs for which A or B level evidence exists:
   - A whole range of drugs have wrong or inconsistent names in VICC and are rectified by SERVE.
   - VICC does not explicitly model the difference between "multiple different drugs" and a "combination treatment of multiple drugs". 
   This gets rectified by SERVE.

### DoCM Curation   

DoCM is used exclusively for known hotspot generation. The filtering is therefore tailored for hotspots:
 - Entries implying general codon mutations are removed.
 - Unusual notation for inframe deletions and insertions is removed.
 - A small number of mutations don't exist on the transcript specified by DoCM and are removed.
 
### iClusion Curation

iClusion contributes to actionability only. SERVE configures every trail to B-level evidence with responsive direction. 
SERVE only considers trials with one or more molecular inclusion criterium, and applies the following filtering and curation:
 - Mutations that are inconsistent with the Hartwig driver catalog are removed
 - Mutations that are detectable in DNA are removed. One example is "Expression" of a gene.
 
Finally, similar to VICC, trials for which no DOIDs have been specified are rectified by SERVE. 
 
## Overview of the SERVE algorithm

SERVE starts with reading the various knowledgebases and after filtering and curation, knowledge is extracted.
A knowledgebase can contribute to known and/or actionable events. Current configuration as follows:
      
Knowledgebase  | Contributes to Known events? | Contributes to Actionable events?
---|---|---
DoCM | Yes | No 
iClusion | No | Yes
VICC | Yes | Yes

Knowledge extraction is performed on a per-knowledgebase level after which all events are consolidated as follows:
 - All known events are aggregated on a per-event level where every level can have a set of sources in which the event has been defined as pathogenic.
 - All actionable events are concatenated so every actionable event that is present in multiple sources will be present multiple times in 
 the actionable output. 




  
 
 
 
  

 
 
 
 
     


  
 
