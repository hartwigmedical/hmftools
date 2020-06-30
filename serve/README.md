# SERVE

**S**earch **E**xternal **R**esources for **V**ariant **E**vidence 

SERVE harmonizes different sources of evidence into a single unified model that can be readily used in genomic analyses.  
The model provides a mapping from genomic events to clinical evidence (responsiveness or resistance for specific drugs). In addition, SERVE
provides lists of genomic events known to be able to drive cancer.  

### Inputs

The following inputs are supported in SERVE:
 - [CGI](https://www.cancergenomeinterpreter.org) - general purpose knowledgebase that is supported through [VICC](http://cancervariants.org)
 - [CIViC](https://civicdb.org) - general purpose knowledgebase that is supported through [VICC](http://cancervariants.org)
 - [CKB](https://ckb.jax.org) - general purpose knowledgebase that is supported through [VICC](http://cancervariants.org)
 - [OncoKB](https://www.oncokb.org) - general purpose knowledgebase that is supported through [VICC](http://cancervariants.org)
 - [iClusion](https://iclusion.org) - a database with all active and recruiting clinical trials in the Netherlands 
 along with cancer types and molecular inclusion criteria for these trials.
 - [CBG Compassionate Use](https://www.cbg-meb.nl/onderwerpen/hv-compassionate-use-programma/overzicht-goedgekeurde-cup) - 
 a database of approved compassionate use programs in the Netherlands
 - HMF Knowledgebase - a database of known driver mutations generated from our database or curated from literature by ourselves.
  
Do note that SERVE does not provide the actual data that is input to the algorithm, but only provides support for its ingestion.
Also, while SERVE itself is open-source, the sources that can be ingested come with their own licenses and up to the user to make sure 
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
 - Single gene events such as amplifications or general inactivation of a gene
 - Types of mutations in genic ranges such as:
    - "inframe insertions in EGFR exon 20"
    - "splice site mutations in MET exon 14"
    - "any type of missense mutation in BRAF codon 600"
 - Specific missense mutations such as BRAF V600E
 
 Missense mutations and genic ranges come in gDNA format. 
 SERVE can be configured to generate its output for either HG19 or HG38 reference genome version
 
 
     


  
 