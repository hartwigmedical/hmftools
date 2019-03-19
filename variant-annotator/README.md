# Variant Annotator

## Introduction

This tool loads Structural Variants from a Structural Variant VCF file into the database. Its primary support is for GRIDSS VCFs, but will temporarily provide backwards-compatibility for Manta VCFs.

The SVs are enriched with Purple copy number data prior to insertion into the structuralVariant table.

It then looks for disruptions and fusions and inserts these into the following tables:
* structuralVariantBreakend - for any SV breakend which falls within a gene (or upstream of the promotor region)
* structuralVariantFusion - for detected fusions, whether reportable or not
* structuralVariantDisruption - for potential disruptions, whether reportable or not

For a full explanation of the fusion-detection logic, see the Hartwig Fusion Detection Logic Google Doc.

Transcript, Exon, Gene and Protein data for fusion logic is sourced from the public Ensembl database, but for efficiency the VariantAnnotator relies on a local cache of the required Ensembl data in the form of 3 files:
* ensembl_gene_data.csv
* ensembl_protein_features.csv
* ensembl_trans_exon_data.csv

These files can be generated in a directory (that needs to exist) by running the following command:

```
java -jar variant-annotator-x.y-with-dependencies.jar
	-write_ensembl_cache \
	-data_output_dir /path/to/ensembl_data_cache \
	-ensembl_db mysql://ensembldb.ensembl.org:3337/homo_sapiens_core_89_37 \
	-ensembl_user anonymous
```

## Resources

|<img src="http://cancer.sanger.ac.uk/images/banners/cosmic_banner_328x68.png" width=148 height=31/>|[Cosmic Gene Fusion CSV](https://www.dropbox.com/s/ettsvttgrg1lc6j/cosmic_gene_fusions.csv?dl=0)|[from: cancer.sanger.ac.uk](http://cancer.sanger.ac.uk/cosmic)|[Curr Protoc Hum Genet 10:.11 (2016)](https://www.ncbi.nlm.nih.gov/pubmed/27727438)|
|---|---|---|---|

## Dependencies

* The MySQL HMF patient-db schema loaded (see patient-db project)
* The Ensembl data cache as described above
* PURPLE data either sourced from the DB or local files

## Usage

```
java -jar variant-annotator-x.y-with-dependencies.jar
    -vcf_file /path/to/bpi.vcf
    -fusion_pairs_csv /path_to/fusion.csv
    -promiscuous_five_csv /path_to/knownPromiscuousFive.csv 
    -promiscuous_three_csv /path_to/knownPromiscuousThree.csv
    -ref_genome /path_to/Homo_sapiens.GRCh37.GATK.illumina.fasta
    -ensembl_data_dir /path/to/ensembl_data_cache/
    -db_url "mysql://localhost/hmfpatients"
    -db_user username
    -db_pass password
    -sample SAMPLE_ID or * to run in batch mode
``` 

Other options for running:
* source_svs_from_db - skip VCF file reading and SV loading, and instead download SVs from the DB and only run fusion and disruption logic
* data_output_dir - write fusion and disruption data to CSV files in this directory 
* skip_db_upload - run fusion and disruption logic but don't upload the results to the DB. This is useful for dry-run testing.
* sample '*' - used in conjunction with source_svs_from_db, will download all samples' SVs from the DB and run fusion & disruption logic for them
* write_ensembl_cache - generate Ensembl data cache files. Need to supply connection details: ensembl_db, ensembl_user and ensembl_pass
* log_debug - verbose logging
