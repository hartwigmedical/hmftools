# **S**earch **E**xternal **R**esources for **V**ariant **E**vidence (SERVE)

SERVE harmonizes various sources of evidence into a single unified model that can be readily used in genomic analyses.  
The model provides a mapping from genomic events to clinical evidence. 
In addition, SERVE generates output containing all genomic events that are implied to be able to driver cancer.  

### Supported input knowledgebases

SERVE supports the ingestion of the following knowledgebases:
 - [CGI](https://www.cancergenomeinterpreter.org) - general purpose knowledgebase that is supported through [VICC](http://cancervariants.org)
 - [CIViC](https://civicdb.org) - general purpose knowledgebase that is supported through [VICC](http://cancervariants.org)
 - [CKB](https://ckb.jax.org) - general purpose knowledgebase that is supported through [VICC](http://cancervariants.org)
 - [OncoKB](https://www.oncokb.org) - general purpose knowledgebase that is supported through [VICC](http://cancervariants.org)
 - [iClusion](https://iclusion.org) - a database with all active and recruiting clinical trials in the Netherlands 
 along with cancer types and molecular inclusion criteria for these trials.
 - [CBG Compassionate Use](https://www.cbg-meb.nl/onderwerpen/hv-compassionate-use-programma/overzicht-goedgekeurde-cup) - 
 a database of approved compassionate use programs in the Netherlands
 - HMF Knowledgebase - a database of known driver mutations generated from the Hartwig database or curated from literature by the Hartwig team.
  
Do note that SERVE does not provide the actual data that is input to the algorithm, but only provides support for its ingestion.
While SERVE itself is open-source, the sources that can be ingested have their own licensing and up to the users to make sure 
they are compliant with the usage of the data itself. 

### Outputs

Clinical evidence is generated in the following datamodel:
 - Drug (along with a drug type)
 - Cancer type for which the evidence applies.
 - Tier / Evidence level
 - Responsive or resistant
 
The following genomic events can be mapped to clinical evidence:
 - Genome-wide events such as signatures or MSI status
 - Multi-gene events such as gene fusions
 - Single gene events such as amplification or general inactivation of a gene
 - Types of mutations in genic ranges such as:
    - inframe insertions in EGFR exon 20
    - splice site mutations in MET exon 14
    - any type of missense mutation in BRAF codon 600
 - Specific missense mutations such as BRAF V600E
 
 Missense mutations and genic ranges come in gDNA format and SERVE can be configured to generate its output for either HG19 or HG38 reference genome version.
 
 In addition to generating a mapping from various genomic events to clinical evidence, 
 SERVE also generates the following outputs describing genomic events implied to be able to driver cancer:
  - Specific known fusion genes including individual genes that are believed to be able to drive cancer when fused with any other gene (in either 5' or 3' position)
  - Known amplifications and deletions
  - Known hotspots (missense and nonsense mutations)

 ## Extraction of genomic events from knowledgebases
 
 ### Protein resolving for missense and nonsense mutations
 
 Evidence on missense and nonsense mutations generally come in their protein annotated form (eg BRAF V600E). SERVE uses [transvar](https://github.com/zwdzwd/transvar) 
 to resolve these annotations into genomic coordinates (referred to as hotspots) for the reference genome version we are generating a model for.
 
 The first step is to choose what transcript to use for converting protein annotation back to genomic coordinates. The rule here is:
  1. If the knowledgebase configured a transcript for a mutation we exclusively use that transcript.
  1. If no transcript is configured we use our internally defined canonical transcript (which is generally the canonical transcript defined by ensembl)
  1. If a protein annotation does not exist on the canonical transcript and has no transcript configured, we use any random transcript for which the protein annotation does exist.
 
 Assuming we found a suitable transcript, we then derive N hotspots for each protein annotation as follows:
  - In case the mutation is caused by SNV or MNV we generate every possible trinucleotide mutation combination that leads to the protein annotation.
    - Exception is in case the codon that includes the SNV or MNV spans multiple exons. For simplicity we only generate a hotspot for the actual SNV/MNV.
  - In case the mutation is caused by a duplication (DUP) or an inframe deletion (DEL) we generate 1 hotspot where we 
  assume the exact reference sequence has been duplicated or deleted. 
    - In case the duplication or deletion spans multiple exons we skip hotspot generation.
  - In case the mutation is caused by an inframe insertion (INS) we have 2 flavours:
    1. In case the insertion is 1 amino acid, we generate hotspots for every trinucleotide coding for that amino acid.
    1. In case the insertion is more than 1 amino acid, we generate just one of the potentially many hotspots.
  - In case of a complex deletion/insertion (DELINS) we extrapolate from the rules for hotspot generation for deletions and insertions. 
  Hence, we assume the reference sequence has been deleted, and one new hotspot is inserted unless the insertion is 1 amino acid in which
  case we generate hotspots for all trinucleotides coding for that amino acid.   
 
 
 
 
     


  
 