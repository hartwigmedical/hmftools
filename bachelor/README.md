# Bachelor

Bachelor filters and annotates all nonsense, splice and frameshift germline variants which affect our included gene list but are not in the blacklist. It also picks up germline variants of other types which are specifically whitelisted.

The steps in the routine are as follows:
1. Parse germline VCF to identify candidate variants. Apply black and white lists to exclude or include variants.
2. If the tumor alt and read depth counts are not in the germline VCF, obtain them from the sample's tumor BAM file.
3. Enrich variants with copy number data
4. Write final set of germline variant data to DB (germlineVariant table) and a TSV output file.

### Whitelist
Any missense inframe indel, synonymous or non-coding variant that is present in Clinvar with a significance IN (‘PATHOGENIC’,’LIKELY PATHOGENIC’) or which has conflicting evidence but none of them ‘BENIGN’ or ‘LIKELY_BENIGN’. 

Any missense variant with a protein coding change that exactly matches the protein coding change of a variant with a significance in clinvar IN (‘PATHOGENIC’,’LIKELY PATHOGENIC’) or which has conflicting evidence but none of them ‘BENIGN’ or ‘LIKELY_BENIGN’. 

### Blacklist
 - Any nonsense, splice or frameshift variant which is present in clinvar with a significance IN (‘BENIGN’,’LIKELY BENIGN’)
 - Gene = BRCA2 + Codon > 3326 
 - CodingEffect in (‘SPLICE’,’NONE’) AND Type = ‘INDEL’ AND (repeatCount>=8 OR indelSequence == microhomology)
 - Any frameshift variant which is offset by another frameshift variant to be net inframe in more than 50% of samples which it is found in our cohort.  These will appear as explicit variants in our blacklist configuration 
 - TSC2 16:2137924 T>TCCCTGCAGTGCAGGAAAGGTAGGGCCGGGTGGGG (rs137854222) - This is classified as a frameshift variant by snpeff but overlaps the splice region and has perfect microhomology and so has no net effect on coding or splicing.

Custom white and blacklistings can be specified per gene in the XML config by either:
- MinCodon OR
- Ref, Alt,Chromosome & Position


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
- include_vcf_filtered - process non-passing variants from the VCF, mark as filter = VCF_FILTERED
- skip_enrichment - 
- log_debug - log in a verbose manner


## Input Files

Bachelor requires that VCF files have been annotated using SnpEff.

## XML Configuration

Bachelor uses 2 config files:
* an XML config file specifying the genes and effects to filter on, and blacklist entries for known exceptions.
* optionally a file of ClinVar data to additionally black and white list variants

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
    -create_filter_file /path/to/clinvar_snp.vcf.gz 
    -output_dir /path/to/dir_that_will_hold_all_output
```
