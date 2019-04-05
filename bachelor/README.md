# Bachelor

Bachelor will report all nonsense, splice or frameshift variants which affect our included gene list but are not in the blacklist.   It will also report variants of other types which are specifically whitelisted.


### Whitelist

Any missense inframe indel, synonymous or non-coding variant that is present in clinvar with a significance IN (‘PATHOGENIC’,’LIKELY PATHOGENIC’)  or which has conflicting evidence but none of them  ‘BENIGN’ or ‘LIKELY_BENIGN’. 

Any missense variant with a protein coding change that exactly matches the protein coding change of a variant with a significance in clinvar IN (‘PATHOGENIC’,’LIKELY PATHOGENIC’) or which has conflicting evidence but none of them ‘BENIGN’ or ‘LIKELY_BENIGN’. 

### Blacklist
 - Any nonsense, splice or frameshift variant which is present in clinvar with a significance IN (‘BENIGN’,’LIKELY BENIGN’)
 - Gene = BRCA2 + Codon > 3326 
 - CodingEffect in (‘SPLICE’,’NONE’) AND Type = ‘INDEL’ AND (repeatCount>=8 OR indelSequence == microhomology)
 - Any frameshift variant which is offset by another frameshift variant to be net inframe in more than 50% of samples which it is found in our cohort.  These will appear as explicit variants in our blacklist configuration 
TSC2 16:2137924 T>TCCCTGCAGTGCAGGAAAGGTAGGGCCGGGTGGGG (rs137854222) - This is classified as a frameshift variant by snpeff but overlaps the splice region and has perfect microhomology and so has no net effect on coding or splicing.

## Usage

Single-sample mode:

```bash
java -jar bachelor.jar 
    -sample [sampleId] 
    -output_dir ./logs 
    -configXml bachelor.xml 
    -ext_filter_file BACHELOR_CLINVAR_FILTERS.csv 
    -runDirectory [path_to_sample_vcf] 
    -log_debug
```

Batch mode:

```bash
java -jar bachelor.jar 
    -sample [sampleId] 
    -output_dir ./logs 
    -configXml bachelor.xml 
    -ext_filter_file BACHELOR_CLINVAR_FILTERS.csv 
    -batchDirectory [path_to_sample_directories] 
    -log_debug
```

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