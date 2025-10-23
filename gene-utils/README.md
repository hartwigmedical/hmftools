# HMF Gene Utilities

This module has routines to generate resources files from Ensembl for use by many of the HMF applications. 

Schema of the Ensembl database can be found on https://m.ensembl.org/info/docs/api/core/core_schema.html


### Generating cached Ensembl data files
The Ensembl data cache is primarily used by Linx, SERVE and Isofox and contains gene, transcript, exon, splice distance and protein data for all genes.

To generate the 4 files, run the following command specifying the Ensembl database instance and the ref-genome version (V37 or V38):

```
java -cp gene-utils.jar com.hartwig.hmftools.geneutils.ensembl.GenerateEnsemblDataCache
    -ensembl_db [ensembl_URL, see below] -ensembl_user [ensembl_user] -ensembl_pass [password]
    -ref_genome_version [37 or 38] 
    -output_dir /path_to_write_data_files/ 
```

The latest Ensembl database URLs for 19/37 & 38 are:
- mysql://ensembldb.ensembl.org:3337/homo_sapiens_core_89_37
- mysql://ensembldb.ensembl.org:3306/homo_sapiens_core_104_38

(connect as 'anonymous' and without a password to access these public Ensembl databases)

Note that ENST00000467125 is blacklisted from Ensembl as it is shares a splice boundary with a chimeric pathogenic GOPC_ROS1 fusion transcript.

### Generating the Ref-Seq data file
Generate the RefSeq gene data file:
```
java -cp gene-utils.jar com.hartwig.hmftools.geneutils.ensembl.GenerateRefSeq
    -ensembl_db [ensembl_URL, see below] -ensembl_user [ensembl_user] -ensembl_pass [password] 
    -output_dir /path_to_write_data_files/ -ref_genome_version [37 or 38]
```


### Gene Mapping between Ref Genomes
A routine to map gene definitions between ref-genome 37 and 38 can be run. It will generate a file mapping GeneId and GeneName between the 2 versions.
Gene matches are tried in the following order:
- GENE_ID - if unchanged
- GENE_NAME - if GeneId ID has changed but GeneName remains the same
- SYNONYM - genes have the same Entrez (typically HGNC) IDs
- COORD - if a lift-over BED file is provided, then a coordinate match is attempted

```
java -cp gene-utils.jar com.hartwig.hmftools.geneutils.mapping.EnsemblGeneMapper
    -ensembl_dir_37 [path_to_ensembl_data_cache_v37] -ensembl_dir_38 [path_to_ensembl_data_cache_v38]
    -lift_over_file [genes_37_lift_over_to_38.bed] 
    -output_dir /path_to_write_mapping_output/
```

## Overview of Gene configuration

### Global Gene Panel
* HMF’s gene model is keyed on HGNC symbol (referred to as “gene” in our data model)
* HMF’s universe of genes consists of all HGNC symbols (http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt) that are matched to an ensembl gene id  that also exists in the latest ensembl version, or failing that gene with matching HGNC symbol in Ensembl.  Notes:
  + If the Ensembl gene name and HGNC gene symbol disagree we use HGNC Symbol as gene name
  + RP11-356O9.1 (hg37)  (annotated as ‘novel transcript in HG38, also known as AL121790.1) is a highly recurrent 5’ fusion partner in Prostate cancer, but is not present in the HGNC so is added to the HGNC symbols and mapped to ENSG00000258414
  + HIST1H3B (hg37) has a different HGNC symbol and Ensembl gene id in HG38 and is enforced to match H3C2
  + The Ensembl gene ids for UGT1A1 and UGT1A8 are swapped between the HG37 and HG38 versions of Ensembl. To resolve this UGT1A1 (hg37) is enforced to match UGT1A1 (hg38), and UGT1A8 (hg37) is enforced to match UGT1A8 (hg38).
  + If multiple Ensembl genes match the same HGNC symbol, then always choose one only: the matching Ensembl gene id first or if none match then the longest transcript.
* HMF’s universe of transcripts consists of all ensembl transcripts belonging to a gene with a matching HGNC symbol
* HMF uses the ensembl canonical transcript as its canonical transcript.   Functional annotations of point mutations are made relative to this transcript, drivers are called on this transcript, and the transcript is favored in fusion prioritisation.  Specific additional transcripts of clinical interest (eg. CDKN2Ap14ARF) may be configured for specific genes in the driver gene panel and will be also used to call drivers in this gene.   Notes:
  + The ensembl method for choosing canonical transcripts has changed significantly for hg38 since Apr 2021, but is frozen for hg37 so the canonical transcripts may differ between versions.
  + If a variant overlaps multiple genes, any driver panel gene is prioritised first, then the gene with the longest canonical transcript (coding bases).

### Curated gene configuration in Hartwig Pipeline
There are 4 additional places in our system where we manually curate lists of genes which are relevant to our pipeline:

* Driver gene panel - Hartwig’s curated panel of genes which is used to define reportable for both OncoAct and the driver catalog.  
* FusionKb - Curated list of gene pairs that are known to form pathogenic gene fusions and which are reportable in OncoAct.   The panel is also used for tiered filtering sensitivity in ESVEE.  
* PEACH - Configuration of genes for pharmacogenomics, currently only DPYD and UGT1A1
* Germline gene notification - Germline reporting rules for OncoAct only.

### Auto-generated interim gene files
A number of gene related interim configuration files are also used by various applications:
* Fusion Hotspot BEDPE - Defines a set of paired regions of extra sensitivity for ESVEE filtering of SV.   Consists of the entire gene on the 5’ partner and the entire gene as well as an additional 10 kb upstream of any 3’ partner for any KNOWN fusions
* HmfRefCDS.RData - representation of canonical transcripts used by PADDLE (dnds)
* DndsDriverLikelihood files - Resource files generated by PADDLE using dndscv for calculating driver likelihood in PURPLE

### Gene configuration upgrade procedures 

* Global gene name or definition changes?
  + Download definitions from HGNC & matching version of ensembl
  + Run GeneUtils GenerateEnsemblDataCache to generate new files for genes with HGNC entries + special overrides/additions
  + Follow PADDLE procedure to regenerate DNDS / driver likelihood files
  + Rerun PAVE on both germline and somatic variants
* Fusion gene name or configuration changes?
  + Review and update FusionKB google sheet & export to TSV 
  + RegenerateFusion bedpe if fusion genes affected (hmftools/gene_panel/GenerateKnownFusionData.R)
  + Rerun ESVEE
* Rerun PURPLE & LINX
* CUPPA refresh (if fusion OR driver genes affected)
  + Regenerate CUPPA reference files 
  + Rerun CUPPA


