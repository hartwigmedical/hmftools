# HMF Gene Utilities

This module ahas routines to generate resources files from Ensembl for use by many of the HMF applications. 

Schema of the Ensembl database can be found on https://m.ensembl.org/info/docs/api/core/core_schema.html


### Generating cached Ensembl data files
The Ensembl data cache is primarily used by Linx and Isofox and contains gene, transcript, exon, splice distance and protein data for all genes.

To generate the 4 files, run the following command specifying the Ensembl database instance and the ref-genome version (V37 or V38):

```
java -cp gene-utils.jar com.hartwig.hmftools.geneutils.ensembl.GenerateEnsemblDataCache
    -ensembl_db [ensembl_URL, see below] -ensembl_user [ensembl_user] -ensembl_pass [password] 
    -output_dir /path_to_write_data_files/ -ref_genome_version [37 or 38]
```

The latest Ensembl database URLs for 19/37 & 38 are:
- mysql://ensembldb.ensembl.org:3337/homo_sapiens_core_89_37
- mysql://ensembldb.ensembl.org:3306/homo_sapiens_core_104_38

(connect as 'anonymous' and without a password to access these public Ensembl databases)

Note that ENST00000467125 is blacklisted from Ensembl as it is shares a splice boundary with a chimeric pathogenic GOPC_ROS1 fusion transcript.


### Generating Gene Panel and Ref-Seq data files
The hmf-common module contains a Ensembl resource file for use in code - all_genes.tsv (V37 and V38 versions). 
It contains gene, transcript and exon data for about 25K genes based on the following criteria:
- canonical transcript
- protein-coding genes or those with an Entrez ID
- a few other specific instances eg non-canonical transcript of CDKN2A p14arf and C11orf95 - a known fusion candidate in children's ependymoma 

```
java -cp gene-utils.jar com.hartwig.hmftools.geneutils.ensembl.GenerateGenePanelRefSeq
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
