# Pave (Prediction and Annotation of Variant Effects)

Pave annotates a somatic variant VCF with gene and transcript coding and protein effects.

For each impacted transcript it will add the following fields:

Field | Description
---|---
Gene | Ensembl gene ID
GeneName | HGNC gene name
Transcript | Ensembl transcript ID
Effects | list of effects separated by '&'
SpliceRegion | true/false - if variant overlaps with the 8 intronic bases or 3 exonic bases around a splice junction 

For any gene with 1 or more impacts, the following summary data is written:

For each impacted transcript it will add the following fields:

Field | Description
---|---
Gene | HGNC gene name
Transcript | Ensembl canonical transcript ID
CanonicalEffect | list of effects separated by '&'
CanonicalCodingEffect | NONE, SPLICE, NONSENSE_OR_FRAMESHIFT, MISSENSE or SYNONYMOUS
SpliceRegion | true/false - if variant overlaps with the 8 intronic bases or 3 exonic bases around a splice junction
HgvsCodingImpact | HGVS coding impact
HgvsProteinImpact | HGVS protein impact 
OtherReportableEffects |Transcript, HGVS Coding, HGVS Protein, Effects, CodingEffect for other reportable transcripts *
WorstCodingEffect | From all transcripts
GenesAffected | Count of genes which the variant overlaps

* if additional reportabled transcripts are configured in drive panel

## Running Pave

### Mandatory Arguments

Argument | Description 
---|---
sample | Sample ID
vcf_file | Input somatic variant VCF
ref_genome | Reference genome fasta file
ensembl_data_dir | Path to Ensembl data cache directory
output_dir | Output directory for VCF and transcript CSV

NOTE: Lilac handles BAMs which have been sliced for the HLA gene regions.

If a sample's tumor BAM is provided in place of the reference BAM, then Lilac will determine the allele solution from it instead.

### Optional Arguments

Argument | Description 
---|---
ref_genome_version | V37 (default) or V38
output_vcf_file | Specify the output VCF filename
only_canonical | Only annotate impacts on canonical transcripts
write_transcript_csv | Write a detailed CSV file for each impacted transcript


```
java -jar pave.jar 
  -sample SAMPLE_ID
  vcf_file /path_to_somatic_vcf_file/
  -ensembl_data_dir /path_to_ensembl_files/
  -ref_genome /path_to_ref_genome_fasta/
  -ref_genome_version [37 or 38] 
  -output_dir /path_to_write_data_files/ 
```

## Overview and algorithm

PAVE predicts the coding impact,protein impact and coding effect of each variant on every overlapping transcript including for up to 1kb upstream.  The following annotations are added for each affected transcript:

* Gene
* TranscriptId
* HGVSCodingImpact
* HGVSProteinImpact
* Effect
* CodingEffect
* SpliceRegion (T/F - any variant that overlaps within 3 exonic or 8 intronic bases of a splice site)
* Population Frequency (gnomAD)

### Effect and Coding Effect
The effects and codingEffects supported by PAVE are the following:

Effect|Coding effect
---|---
• upstream_gene_variant (<1kb)<br />• intron_variant<br />• 5_prime_UTR_variant<br />• 3_prime_UTR_variant<br />• non_coding_transcript_exon_variant | NONE
• synonymous_variant | SYNONYMOUS
• missense_variant<br />• inframe_insertion<br />• inframe_deletion<br />• phased_inframe_insertion<br />• phased_inframe_deletion | MISSENSE
• stop_gained<br />• frameshift<br />• start_lost<br />• stop_lost | NONSENSE_OR_FRAMESHIFT
• splice_donor_variant (D-1,D+1,D+2,D+5)<br />• splice_acceptor_variant (A+1;A+2; A+3 if ALT=G only) | SPLICE

### HGVS Coding Impact

For coding transcripts use c. and for non coding transcripts use n.

Nucleotide numbering conventions:
* No nucleotide 0.  
* Nucleotide 1 is the A of the ATG-translation initiation codon (For non coding 1 is the 1st base of the transcript)
* Nucleotide -1 is the 1st base prior to the ATG-translation initiation codon
* Nucleotide *1 is the 1st base 3’ of the translation stop codon
* Use most 3’ position in case of homology for duplications and deletions (with an exception for homology at splice boundary.  See https://www.hgvs.org/mutnomen/recs-DNA.html#except)

Examples:



### HGVS Protein Impact

Phased inframe variants should get the combined impact of both variants

Amino acid numbering conventions
* AA 1 = Translation initiator Methionine
* AA *110+1 = 1st AA after the stop codon at AA position 110
* Use most 3’ residue in case of AA homology

Examples:



### Population Frequency

We annotate the population frequency using Gnomad (v3.1.2 for hg38, v2.1.1 for GRCH37).  We filter the Gnomad file for variants with at least 0.00005 frequency and and we annotate with a resolution of 0.0001. 

