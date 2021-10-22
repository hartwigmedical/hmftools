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

