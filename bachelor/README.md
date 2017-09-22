# Bachelor

Bachelor is a tool to check whether a patient (represented by a set of VCF files) is eligible for any number of configured studies / programs.

Programs are found in the configDirectory, and are in XML format. 

## Usage
```bash
java -jar bachelor-1-jar-with-dependencies.jar
    -configDirectory /path/to/xmls
    -vcf germline.vcf
    -vcf somatic.vcf
    -vcf structuralVariants.vcf
```

## Input Files

Bachelor requires that VCF files have been annotated using SnpEff.

## XML Configuration

Bachelor loads all the xml files in the given *configDirectory*.

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