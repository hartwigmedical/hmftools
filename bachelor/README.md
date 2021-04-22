# Note

Bachelor functionality is now deprecated and replaced by relevant functionality in Sage (for directly calling germline variants) and PURPLE (for relevant annotation steps).

# Bachelor

Bachelor evaluates the pathogenicity of germline variants from a snpeff annotated germline vcf in a configured panel of genes. Pathogenicity is determined based both on ClinVar and a list of configured snpeff effects that are deemed to be likely pathogenic even if unannotated in ClinVar.

The steps in the routine are as follows:
1. Parse germline VCF to identify candidate variants. 
2. If the tumor alt and read depth counts are not in the germline VCF, obtain them from the sample's tumor BAM file.
3. Enrich variants with copy number data from purple
4. Enrich variants with ClinVar pathogenicity information and configured whitelist / blacklist and determine reportability
5. Write final set of germline variant data to DB (germlineVariant table) and a TSV output file.

## Annotation
In addition to standard PURPLE annotations (https://github.com/hartwigmedical/hmftools/tree/master/purple#10-somatic-enrichment), Bachelor also annotates the following information:

### Pathogenicity
Pathogenicity is determined primarily from ClinVar, or may be set or overridden via a configured blacklist or whitelist.  

Permitted values are:

* BLACK_LIST - Variant matches black list configuration of known benign variant. INDELs in repeats or with microhomology in splice and intronic regions are also blacklisted (ie. CodingEffect in (‘SPLICE’,’NONE’) AND Type = ‘INDEL’ AND (repeatCount>=8 OR indelSequence == microhomology))
* WHITE_LIST - Variant matches white list configuration of known pathogenic variants
* CLINVAR_PATHOGENIC - At least 1 interpretation of 'PATHOGENIC' and none ‘BENIGN’ or ‘LIKELY_BENIGN’
* CLINVAR_LIKELY_PATHOGENIC - No interpretation of PATHOGENIC, but at least 1 interpretation of 'LIKELY_PATHOGENIC' and none ‘BENIGN’ or ‘LIKELY_BENIGN’
* CLINVAR_CONFLICTING - Variant has both likely 'BENIGN'/'LIKELY_BENIGN' and 'PATHOGENIC'/'LIKELY_PATHOGENIC' interpretations
* CLINVAR_LIKELY_BENIGN - No interpretation of 'BENIGN' and at least 1 interpretation of 'LIKELY_BENIGN' and none ‘PATHOGENIC’ or ‘LIKELY_PATHOGENIC’
* CLINVAR_BENIGN - At least 1 interpretation of 'BENIGN' and none ‘PATHOGENIC’ or ‘LIKELY_PATHOGENIC’
* UNANNOTATED - Variant is not annotated in Clinvar

Note that custom white and blacklists can be specified per gene in the XML config by either:
- MinCodon OR
- Ref, Alt,Chromosome & Position

### Filter
A filter field is populated. Permitted values are:
* PASS
* GERMLINE_FILTERED - Variant was soft filtered in input germline vcf file
* ARTEFACT - Phred Score < 150 and adjusted VAF < 0

### Reported
This is a summary boolean flag which implements the HMF logic for reportable germline variants. A variant will be marked as reported if the filter = PASS and either one or both of the following criteria are met:
* Pathogenicity in ('WHITE_LIST','CLINVAR_PATHOGENIC','CLINVAR_LIKELY_PATHOGENIC')
* Pathogenicity = 'UNANNOTATED' and effect is configured as a known snpeffect

## Usage
1. Find candidate germline variants, write to file and optionally write to db.

```bash
java -jar bachelor.jar 
    -sample [sampleId] 
    -germline_vcf /sample_path/sampleid.germline_variants.vcf.gz
    -tumor_bam_file /sample_path/sampleid.tumor.bam 
    -purple_data_dir /sample_path/purple/
    -xml_config /path/bachelor_config.xml 
    -ext_filter_file /path/bachelor_clinvar_filters.csv
    -ref_genome /path/ref_genome.fasta 
    -output_dir /sample_path/bachelor/ 
    [-db_url [db_url] -db_user [user] -db_pass [password]] 
```

2. Load germline variants from file and upload to database

```bash
java -cp bachelor.jar com.hartwig.hmftools.bachelor.LoadGermlineVariants 
    -sample [sampleId] 
    -sample_data_dir /sample_path/bachelor/
    -db_url [db_url] -db_user [user] -db_pass [password] 
```

Optional config:
- include_vcf_filtered - process non-passing variants from the VCF, mark as filter = GERMLINE_FILTERED
- log_debug - verbose logging


## Input Files
Bachelor requires that VCF files have been annotated using SnpEff.

## Configuration
Bachelor uses 2 config files:
* an XML config file specifying the genes and effects to filter on, and blacklist entries for known exceptions
* optionally a file of ClinVar variant data to set pathogenic status

This very basic example has a single criteria: missense variants in BRCA2.

```xml
<?xml version="1.0" ?>
<Program name="ExampleStudy" xmlns="http://www.hartwigmedicalfoundation.nl/bachelor.xsd">
    <Panel>
        <Gene name="BRCA2" ensembl="ENST00000544455" refseq="NM_000059"/>
        <SnpEffect>missense_variant</SnpEffect>
    </Panel>
</Program>
```


## ClinVar Filter Creation
Create a set of files from ClinVar variant classification with this command:

```bash
java -cp bachelor.jar com.hartwig.hmftools.bachelor.ExternalDBFilters
    -xml_config /path/to/bachelor_config.xml 
    -clinvar_db_filter_file /path/to/clinvar_snp.vcf.gz 
    -output_filter_file /path/to/output_filter_file
```

log_debug can be passed an optional additional parameter for additional logging.
