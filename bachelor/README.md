# Bachelor

Bachelor will report all nonsense, splice or frameshift variants which affect our included gene list but are not in the blacklist.   It will also report variants of other types which are specifically whitelisted.

The steps in the routine are as follows:
1. Parse germline VCF to identify candidate variants. Apply black and white lists to exclude or include variants.
2. If the tumor alt and read depth counts are not in the germline VCF, obtain them from the sample's tumor BAM file.
3. Enrich variants with copy number data
4. Write final set of germline variant data to DB (germlineVariant table) and an CSV output file.

In the case where step 2 needs to be completed separately, Bachelor can be run in 2 stages:
- Stage 1 - step 1 above
- Stage 2 - steps 3 & 4 above 

### Whitelist

Any missense inframe indel, synonymous or non-coding variant that is present in clinvar with a significance IN (‘PATHOGENIC’,’LIKELY PATHOGENIC’)  or which has conflicting evidence but none of them  ‘BENIGN’ or ‘LIKELY_BENIGN’. 

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

The default mode for Bachelor is to run in a continuous mode where all 4 steps are run sequentially.

Alternatively Stage 1 and Stage 2 can be run independently as shown below.

1. Single-sample default/continuous mode:

```bash
java -jar bachelor.jar 
    -sample [sampleId] 
    -xml_config bachelor.xml 
    -ext_filter_file BACHELOR_CLINVAR_FILTERS.csv 
    -sample_data_dir [path_to_sample_vcf] 
    -bam_direct 
    -purple_data_dir purple
    -db_url [db_url] -db_user [user] -db_pass [password] 
    -high_confidence_bed [path_to_file] 
    -ref_genome [path_to_file] 
    -write_to_db  
```

2. Single-sample Stage 1 only

```bash
java -jar bachelor.jar 
    -sample [sampleId] 
    -xml_config bachelor.xml 
    -ext_filter_file BACHELOR_CLINVAR_FILTERS.csv 
    -sample_data_dir [path_to_sample_vcf] 
```

3. Single-sample Stage 2 only
```bash
java -jar bachelor.jar 
    -sample [sampleId] 
    -sample_data_dir [path_to_sample_vcf] 
    -bam_direct 
    -purple_data_dir purple
    -db_url [db_url] -db_user [user] -db_pass [password] 
    -high_confidence_bed [path_to_file] 
    -ref_genome [path_to_file] 
    -write_to_db  
```


4. Batch mode:

```bash
java -jar bachelor.jar 
    -sample "*" or omit 
    -xml_config bachelor.xml 
    -ext_filter_file BACHELOR_CLINVAR_FILTERS.csv 
    -sample_data_dir [path_to_sample_vcf] 
    -bam_direct 
    -purple_data_dir purple
    -db_url [db_url] -db_user [user] -db_pass [password] 
    -high_confidence_bed [path_to_file] 
    -ref_genome [path_to_file] 
    -write_to_db  
```

Optional config:
- log_debug - log in a verbose manner
- batch_output_dir - only applicable if running in batch mode, writes variants to a single CSV output file
- purple_data_dir - if omitted, will retrieve copy number data from DB 



## Input Files

Bachelor requires that VCF files have been annotated using SnpEff.

## XML Configuration

Bachelor uses 2 config files:
* an XML config file specifying the genes and effects to filter on, and blacklist entries for known exceptions.
* optionally a file of Clinvar data to additionally black and white list variants

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


## Clinvar Filter Creation

Create a set of files from CLINVAR variant classification with this command:

```bash
java -cp bachelor.jar com.hartwig.hmftools.bachelor.ExternalDBFilters
    -xml_config bachelor.xml 
    -create_filter_file /path_to/clinvar_snp.vcf.gz 
    -output_dir path_to_log
```
