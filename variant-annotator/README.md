# Variant Annotator

## Introduction

This tool will load Structural Variants from a Manta SV VCF which has been post-processed by BPI.

It will preform annotation of each end of a SV, with gene, transcript and exon information.
It will then provide further annotations of *Gene Disruptions* and *Gene Fusions*.

**Gene Disruptions** are where a single end of a SV disrupts the coding of a gene.

**Gene Fusions** are where both ends of a SV are:
* in introns
* the exons being joined are in-phase

## Resources

|<img src="http://cancer.sanger.ac.uk/images/banners/cosmic_banner_328x68.png" width=148 height=31/>|[Cosmic Gene Fusion CSV](https://www.dropbox.com/s/ettsvttgrg1lc6j/cosmic_gene_fusions.csv?dl=0)|[from: cancer.sanger.ac.uk](http://cancer.sanger.ac.uk/cosmic)|[Curr Protoc Hum Genet 10:.11 (2016)](https://www.ncbi.nlm.nih.gov/pubmed/27727438)|
|---|---|---|---|

## Dependencies

* A MySQL (or possibly MariaDB) which has the HMF patient-db schema loaded (see patient-db project)
* Access to the exact ensembl db version used by hmftools ensembl-db submodule
* For Purity-Adjusted VAFs to be recorded, PURPLE must have been run and persisted to DB prior to this tool.

## Usage

```
java -jar variant-annotator-x.y-with-dependencies.jar
    -ensembl_db "mysql://ensembldb.ensembl.org:3337/homo_sapiens_core_89_37?user=anonymous"
    -fusion_csv /path/to/fusion.csv
    -db_url "mysql://localhost/hmfpatients"
    -db_user username
    -db_pass password
    -sample SAMPLE_ID
    -vcf_file /path/to/bpi.vcf
``` 
