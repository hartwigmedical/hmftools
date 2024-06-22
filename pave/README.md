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
HgvsCodingImpact | HGVS coding impact
HgvsProteinImpact | HGVS protein impact 

For any variant with one or more impacted transcripts, the following summary data is written:

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

<nowiki>*</nowiki> if additional reportable transcripts are configured in driver gene panel

## Running Pave

### Mandatory Arguments

Argument | Description 
---|---
sample | Sample ID
vcf_file | Input variant VCF
ref_genome | Reference genome fasta file
ensembl_data_dir | Path to Ensembl data cache directory
driver_gene_panel|Driver Gene Panel
ref_genome_version | 37 (default) or 38

### Optional Arguments

Argument | Description 
---|---
output_dir | Output directory for VCF and transcript CSV, will use input VCF directory if not specified
output_vcf_file | Specify the output VCF filename
only_canonical | Only annotate impacts on canonical transcripts
read_pass_only | Only process passing variants
threads | Splits variants by chromosome across threads
write_pass_only | Only write passing variants
write_transcript_data | Write a detailed TSV file for each impacted transcript

```
java -jar pave.jar 
  -sample SAMPLE_ID
  -vcf_file /path_to_somatic_vcf_file/
  -ensembl_data_dir /path_to_ensembl_files/
  -driver_gene_panel /path_to_gene_panel/
  -ref_genome /path_to_ref_genome_fasta/
  -ref_genome_version [37 or 38] 
  -output_dir /path_to_write_data_files/ 
  -threads 8
```

### Optional Annotations
The following annotations can be applied by PAVE:

Argument | Description | VCF Tag
---|---|---
pon_file | PON file to annotate variants (see file format below) | PON_COUNT, PON_MAX
pon_filters | Apply PON filters (see details below) | Filter 'PON'
mappability_bed | BED file with mappability values | MAPPABILITY=value
clinvar_vcf | VCF from Clinvar database | Writes CLNSIG=significance and CLNSIGCONF=conflicting info
gnomad_freq_file | CSV with Gnomad frequencies per variant
gnomad_freq_dir | Path to Gnmoad frequency files per chromosome
blacklist_bed | BED file with blacklist entries | BLACKLIST_BED
blacklist_vcf | VCF file with blacklist entries | BLACKLIST_VCF


## Overview and algorithm

PAVE predicts the coding impact, protein impact and coding effect of each variant on every overlapping transcript including for up to 1kb upstream.  Note that PAVE takes into account nearby mutations in the same local phase set (LPS) where they alter effect (either same codon, or multiple frameshifts leading to a single inframe indel).  If a variant exists in multiple LPS then the one with the highest support is used.

The following annotations are added for each affected transcript:

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

Effect<sup>1</sup>|Coding effect<sup>2</sup>
---|---
• upstream_gene_variant (<1kb)<br />• intron_variant<br />• 5_prime_UTR_variant<br />• 3_prime_UTR_variant<br />• non_coding_transcript_exon_variant | NONE
• synonymous_variant<br />• phased_synonymous<sup>4</sup> | SYNONYMOUS
• missense_variant<br />• inframe_insertion<sup>3</sup><br />• inframe_deletion<sup>3</sup><br />• phased_inframe_insertion<sup>4</sup><br />• phased_inframe_deletion<sup>4</sup><br />• phased_missense<sup>4</sup> | MISSENSE
• stop_gained<br />• frameshift<br />• start_lost<sup>5</sup><br />• stop_lost<sup>5</sup> | NONSENSE_OR_FRAMESHIFT
• splice_donor_variant (D-1,D+1,D+2,D+5)<br />• splice_acceptor_variant (A+1;A+2; A+3 if ALT=G only) | SPLICE<sup>6,7</sup>

Notes:
1. Multiple effects are possible and should both be annotated in the following circumstances:
   - An Inframe_insertion, inframe_deletion or frameshift variant may also be stop_gained or stop_lost
   - Any deletion or MNV unambiguously overlapping an exon boundary OR an SNV at splice_donor(D-1) will be both splice_acceptor/spice_donor AND one of synonymous_variant, stop_gained, missense_variant,  5_prime_UTR_variant, 3_prime_UTR_variant or non_coding_transcript_exon_variant
2. Where a variant has multiple effects, PAVE ranks in the following order for codingEffect: 
   - NONSENSE_OR_FRAMESHIFT
   - SPLICE
   - MISSENSE
   - SYNONYMOUS
   - NONE
3. Inframe INDELs may occasionally be annotated as notionally partially or completely outside the coding region due to left alignment and microhomology. Any INDEL with a length that is a multiplier of 3, that can be right aligned to be fully inside the coding regions should be marked as effect=inframe_insertion/inframe_deletion (notable examples include known pathogenic variants in KIT (4:55593579 CAGAAACCCATGTATGAAGTACAGTGGA > C) and EGFR (7:55248980 C > CTCCAGGAAGCCT)). 
4. Where there are 2 or more frameshift variants with the same LPS (local phase set), if the combined impact causes an inframe indel, then mark both as effect = phased_inframe_deletion / phased_inframe_insertion.   If a phased inframe indel and snv or mnv affect the same codon, then mark both as phased_inframe_deletion / phased_inframe_insertion and calculate the combined coding effect (eg.  EGFR p.Glu746_Ser752delinsVal). If the net impact is no inserted or deleted bases, then classify the phased variants as phased missense or synonymous. 
5. Where an INDEL also leads to a stop_lost or start_lost, the lost effects are prioritised
6. A SPLICE MNV needs to be marked as splice if any base overlaps a splice site.
7. Any INDEL which overlaps a canonical splice region (ie.[D-1:D+5] OR [A+3:A+1]) should be marked as splice_donor/splice_acceptor if and only if the canonical sites are changed according to the SPLICE rules listed above. Where an INDEL has microhomology extends over a splice donor or splice acceptor region, the variant is tested at both the leftmost and rightmost alignment, with intronic only effects prioritised highest, then exonic effects and finally splice effects.   A notable recurrent example where D+5 is not affected by an indel with microhomology in GRCH37 are indels at the homopolymer at MSH2 2:47641559.  Both splice and frameshift/inframe effects may be reported together if a deletion unambiguously both overlaps coding bases and changes canonical splice sites.


### HGVS Coding Impact

For coding transcripts use c. and for non coding transcripts use n.

Nucleotide numbering conventions:
* No nucleotide 0.  
* Nucleotide 1 is the A of the ATG-translation initiation codon (For non coding 1 is the 1st base of the transcript)
* Nucleotide -1 is the 1st base prior to the ATG-translation initiation codon
* Nucleotide *1 is the 1st base 3’ of the translation stop codon
* Use most 3’ position in case of homology for duplications and deletions (with an exception for homology at splice boundary.  See https://www.hgvs.org/mutnomen/recs-DNA.html#except)

TERT only is also annotated into the upstream region for 300 bases upstream of the start codon of the canonical transcript due to the widespread conventions around TERT promoter mutations

Examples:

Type | Location | Examples | Notes
---|---|---|---
Substitutions (SNV)  | Coding | c.76A>C  | Includes start and stop codons
_ | 5’UTR | c.-14G>C
_ | Intronic (post) | c.88+1G>T
_ | Intronic (pre) | c.89-2A>C
_ | 3’UTR | c.*46T>A
Deletions | Coding | c.76_78delACT
_ | Intronic (pre) | c.726-5537_726-5536delTT
_ | 5’UTR Intronic | c.-147-1093delA
Duplications | Coding | c.128dupA | Use duplications in case of INS with full homology match.
_ | Intronic (post) | c.830+11459_830+11461dupGGA
_ | Start codon overlap | c.-1_1dupAA
Insertions | Coding | c.1033_1034insA
_ | Intronic (post) | c.15+1619_15+1620insTTTGTT
_ | Intronic (5’UTR) | c.-23-304_-23-303insA
Substitutions (MNV) | Coding | c.1390_1391delCGinsTT
_ | Intronic | c.853-2260_853-2258delCACinsTAT
Complex Indels | Coding | c.112_117delinsTG

### HGVS Protein Impact

Phased inframe variants should get the combined impact of both variants

Amino acid numbering conventions
* AA 1 = Translation initiator Methionine
* AA *110+1 = 1st AA after the stop codon at AA position 110
* Use most 3’ residue in case of AA homology

Any variant with codingEffect = SPLICE in a coding transcript is given the hvgs protein impact 'p.?'.  Other variants are annotated as per below examples:

Type | Context | Examples | Notes
---|---|---|---
SYNONYMOUS | Synonymous | p.Leu54= | Snpeff uses p.Leu54Leu but the recommendation has changed (https://www.hgvs.org/mutnomen/disc.html#silent)
_ | Synonymous (MNV multiple codon) | p.Leu54_Arg55=
_ | Synonymous (Stop retained) | p.Ter26= |
MISSENSE | Missense | p.Trp26Cys
_ | Missense (MNV multiple codon) | p.Ala100_Val101delinsArgTrp  | SnpEff uses format: p.AlaVal100ArgTrp
NONSENSE OR FRAMESHIFT  | Stop Gained | p.Trp26*
_ | Stop Gained (MNV multiple codon - 2nd codon stop) | p.Cys495_Val496delinsArg*
_ | Stop Gained (MNV multiple codon - 1st codon stop) | p.Cys495_Val496delins*
_ | Frameshift | p.Arg97fs | Always use simply fs even if also stop gained
_ | Stop Lost | p.Ter407Trpext*? | Ie. a STOP at AA407 changes to a W and extends the protein
_ | Start Lost | p.Met1? | Any variant that disrupts initiator codon
INFRAME | Deletion (single AA) | p.Lys2del
_ | Deletion (range) | p.Gly4_Gln6del 
_ | Deletion (non conservative) | p.Cys28_Lys29delinsTrp
_ | Duplication (single AA) | p.Gln8dup
_ | Duplication (range) | p.Gly4_Gln6dup
_ | Insertion | p.Lys2_Leu3insGlnSer
_ | Insertion (conservative stop) | p.Ser81_Val82ins* | Ie. a STOP codon is inserted between AA81 and AA82
_ | Insertion (non conservative) | p.Cys28delinsTrpVal 
MIXED | Inframe Deletion with stop lost (single AA) | p.104Terdelext*?
_ | Inframe Deletion with stop lost (multiple AA non conservative) | p.Val98_Ter104delinsArgext*?
_ | Inframe insertion + stop gained | p.Leu339delinsHisPhe* | Ignore any AA inserted AFTER the stop codon

### PON Annotation and Filtering
Pave can annotate with PON values if the config 'pon_file' is used. The PON file must be a TSV with the following fields:

```
Chromosome      Position        Ref     Alt     SamplesCount    MaxSampleReads  TotalReads
1       10003   A       T       10      30      191
1       10006   C       A       16      11      86
1       10007   T       G       40      16      234
```

Pave will then add VCF tags 'PON_COUNT' from SamplesCount and 'PON_MAX' from MaxSampleReads for any matched variant.

If the config 'pon_filters' is used, then Pave will additionally add the filter 'PON' to any variant which exceeds both the specified SamplesCount and MaxSampleReads values. 
The filters can be set per variant tier in the form: 'TIER;SAMPLE_COUNT_LIMIT;MAX_READS_LIMIT, for example

```HOTSPOT:5:5;PANEL:2:5;UNKNOWN:2:0```

will mark any variant of tier = HOTSPOT as PON if it matches an entry with 5+ SamplesCount and 5+ MaxSampleReads, 2+ and 5+ for a PANEL tier variant, and 2+ for any other tier variant.

### GNOMAD Population Frequency
We annotate the population frequency using gnomAD v3.1.2 for hg38 (merged with gnomAD v2.1.1 liftover for exome regions only) and v2.1.1 exome only for GRCH37. We filter the Gnomad file for variants with at least 1e-5 frequency for exome only and 5e-5 for genome. The VCF tag 'GND_FREQ' will report the frequency.

### CLINVAR
If a clinvar VCF is provided, PAVE also annotates the clinical signficance of each variant.

### PON settings used in the HMF pipeline

A summary of the PON annotation and filtering currently used in the HMF pipeline is below:

Filter | Annotations | Source | Filter Thresholds | Ref Genome versions
---|---|---|---|---
PON_GNOMAD | GND_FREQ | Gnomad v3 | GND_FREQ<0.00015 | 38 only
PON_PANEL_ARTEFACT | PON_PANEL | Curated FFPE Panel Artefacts*** | PON_PANEL {ANY} | 38 only 
PON | PON_COUNT* | PON_MAX** | HMF Cohort | See detailed table below | 37 & 38

<nowiki>*</nowiki> Count germline samples with at least 3 reads and sum of base quality > 30
** Maximum read support in any one sample
*** The FFPE panel artefacts were curated from recurrent variants in the panel regions only of 35 samples run on HMF FFPE tumor only panels.

The filters for the HMF cohort PON depend on the ref genome version and are as follows:

Version | Samples | HOTSPOT | PANEL | OTHER
---|---|---|---|---
37 | 1000 | PON_MAX>=5 & PON_COUNT>=10 | PON_MAX>=5 & PON_COUNT>=6 | PON_COUNT>=6
38 | 98 | PON_MAX>=5 & PON_COUNT>=5 | PON_MAX>=5 & PON_COUNT>=2 | PON_COUNT>=2

## Known Issues
- Frameshifts may not always be fully aligned to 3'UTR for HGNC protein annotation
- Where multiple ALTs are included on a single line only the 1st ALT allele will be annotated.   A workaround is to split multiallelic lines into multiple records first (eg.  the -m none option in bcftools).
- Duplications may sometimes be marked as insertions 
- Inserts between the translation/coding start base and the previous base may not be classified correctly and/or have the wrong codon impact due to left alignment

# Version History and Download Links
- [1.4](https://github.com/hartwigmedical/hmftools/releases/tag/pave-v1.4.2)
- [1.3](https://github.com/hartwigmedical/hmftools/releases/tag/pave-v1.3)
- [1.2](https://github.com/hartwigmedical/hmftools/releases/tag/pave-v1.2)
- [1.1](https://github.com/hartwigmedical/hmftools/releases/tag/pave-v1.1)
- [1.0](https://github.com/hartwigmedical/hmftools/releases/tag/pave-v1.0)
