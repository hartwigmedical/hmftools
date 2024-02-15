## Running SAGE in Germline Mode

SAGE can be configured to detect germline variants by switching the tumor and germline samples, and adjusting filtering parameters, panel definition and hotspot inputs.  

If a 'coverage_bed' file is provided, then calculate and write coverage statistics per gene in the PANEL file. This can be used to estimate the likelihood of missing a genuine germline variant.

The following sections describe how SAGE can be configured to detect pathogenic germline variants. 

## Panel
The panel is constructed from the supplied DriverGenePanel.xx.tsv file using [Gene Utils](../gene-utils/README.md). 

All germline and somatic reportable genes are included in the germline panel. 
Unlike the somatic panel, the germline panel includes the UTR regions in addition to the coding regions. 
Splice sites (+1,+2,+5,-2,-1) are included as well.
 
## Hotspots / Whitelist
The germline hotspot is created at the same time as the germline panel. 
We select all variants from the clinvar file that are Pathogenic or Likely_pathogenic (but not Benign or Likely_benign) and are marked as report germline hotspot in the driver gene panel.

A select number of [whitelist](../gene-utils/src/main/resources/drivers/GermlineHotspots.whitelist.38.vcf) variants are added to the hotspot file. 

## BlackList
[Blacklisted](../gene-utils/src/main/resources/drivers/GermlineHotspots.blacklist.38.vcf) variants are not reported as pathogenic.  
Neither are BRCA2 variants in the range 13:32972625-32972907 (or v38 -> chr13:32398488-32398770).

## Parameters
To run SAGE in germline mode we use the following parameters:

```
-tumor REFERENCE_SAMPLE
-tumor_bam REFERENCE_BAM
-reference TUMOR_SAMPLE
-reference_bam TUMOR_BAM
-hotspot_min_tumor_qual 50
-panel_min_tumor_qual 75
-ref_sample_count 0
-panel_only
-panel_bed /opt/resources/sage/37/ActionableCodingPanel.37.bed.gz
-coverage_bed /opt/resources/sage/37/CoverageCodingPanel.37.bed.gz 
-high_confidence_bed /opt/resources/giab_high_conf/37/HighConfidence.37.bed.gz 
-hotspots /opt/resources/sage/37/KnownHotspots.germline.37.vcf.gz 
-ref_genome_version 37 
-ref_genome /opt/resources/reference_genome/37/Homo_sapiens.GRCh37.GATK.illumina.fasta
-ensembl_data_dir /path_to_ensmebl_cache/ \
-out /data/output/TUMOR_SAMPLE.sage.germline.vcf.gz 
``` 

Note that the tumor and reference labels/bams are switched. 

These changes (specifically ref_sample_count=0) disable the germline filters (which is actually the tumor).

## Post Process

Run Pave (https://github.com/hartwigmedical/hmftools/tree/master/pave) with germline reference data and settings:

```
java -jar pave.jar 
  -sample SAMPLE_ID
  -vcf_file /path_to_germline_vcf/
  -ensembl_data_dir /path_to_ensembl_files/
  -driver_gene_panel /ref_files/DriverGenePanel.37.tsv \
  -ref_genome /path_to_ref_genome_fasta/
  -ref_genome_version [V37 or V38] 
  -mappability_bed /ref_files/mappability_150.37.bed.gz \
  -clinvar_vcf /ref_files/clinvar.37.vcf.gz \
  -blacklist_bed /ref_files/KnownBlacklist.germline.37.bed \
  -blacklist_vcf /ref_files/KnownBlacklist.germline.37.vcf.gz \
  -read_pass_only \
  -output_dir /output_path/ 
```

## Known Issues

INDEL VAFs in long repeats are systematically underestimated by SAGE (due to jitter).
This makes things even worse when evaluating variant copy number in tumors for germline variants as we will assume germline VAF = 0.5 if we make the germline genotype HET incorrectly.
