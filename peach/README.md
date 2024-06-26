# PEACH

**P**harmacogenomic **E**valuator **A**nd **C**aller of **H**aplotypes (PEACH) is a pharmacogenomics tool developed for
the [Hartwig Medical Foundation pipeline](https://github.com/hartwigmedical/pipeline5).
It imports haplotypes and related variants from config files, reports the presence of these variants in a
germline VCF, and infers the simplest combination of haplotypes that explains the presence of these variants.

An earlier Python version of this tool is available [here](https://github.com/hartwigmedical/peach).

## Contents
* [Installation](#installation)
* [Arguments](#arguments)
  + [Example Usage](#example-usage)
  + [Mandatory Arguments](#mandatory-arguments)
  + [Optional Arguments](#optional-arguments)
* [Input](#input)
  + [VCF](#vcf)
  + [Haplotype TSV](#haplotype-tsv)
  + [Function TSV](#function-tsv)
  + [Drug TSV](#drug-tsv)
* [Output](#output)
  + [Events TSV](#events-tsv)
  + [Gene Events TSV](#gene-events-tsv)
  + [All Haplotypes TSV](#all-haplotypes-tsv)
  + [Best Haplotypes TSV](#best-haplotypes-tsv)
  + [GQC TSV](#qc-tsv)
* [Algorithm](#algorithm)
  + [Interpret VCF Variant Calls](#interpret-vcf-variant-calls)
  + [Determine Valid Haplotype Combinations](#determine-valid-haplotype-combinations)
  + [Determine the Best Haplotype Combination](#determine-the-best-haplotype-combination)
  + [each QC status](#peach-qc-status)
  + [Examples](#examples)
* [Known Issues / Points for Improvement](#known-issues--points-for-improvement)
* [Version History / Download Links](#version-history--download-links)

## Installation
To install, download the latest compiled jar file from the [download links](#version-history--download-links).

PEACH requires Java 11+.

TODO: Add link to config files

## Arguments
### Example Usage
```
java -jar "/path/to/peach.jar" \
  -vcf_file "/path/to/input.vcf.gz" \
  -sample_name "SAMPLE1" \
  -haplotypes_file "/path/to/haplotypes.tsv" \
  -output_dir "/path/to/output_directory"
```

### Mandatory Arguments
| Long&nbsp;Argument | Description                                                                            |
|--------------------|----------------------------------------------------------------------------------------|
| --vcf_file         | Path to germline VCF file of sample. For instance the germline VCF output from PURPLE. |
| --haplotypes_file  | Path to the TSV that configures the haplotypes that can be called.                     |
| --sample_name      | Sample name. This value should correspond to the ID in the input VCF file.             |
| --output_dir       | Directory to write the output to.                                                      |

### Optional Arguments

| Long&nbsp;Argument | Description                                                                                                                                          |
|--------------------|------------------------------------------------------------------------------------------------------------------------------------------------------|
| --function_file    | Path to TSV describing function of configured haplotypes.                                                                                            |
| --drugs_file       | Path to TSV with relevant drugs for configured genes and URLs to suggested prescription adjustments for these drugs depending on haplotype function. |

## Input
### VCF
PEACH has been designed to work with VCF files that follow the VCF Version 4.2 format, see
[specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf). 
One of the samples in the VCF should have a label equal to the configured sample name,
and the "GT" subfield for this sample should be included and filled in with diploid calls.

### Haplotype TSV
| Column    | Example                                       | Description                                                                                                                                                                                                                                   |
|-----------|-----------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| gene      | `DPYD`                                        | Name of gene.                                                                                                                                                                                                                                 |
| haplotype | `*B3`                                         | Name of haplotype.                                                                                                                                                                                                                            |
| default   | `false`                                       | `true` if this is the haplotype of this gene in the reference genome. `false` otherwise.                                                                                                                                                      |
| wildType  | `false`                                       | `true` if this haplotype is considered to be the wild-type haplotype of this gene. `false` otherwise.                                                                                                                                         |
| events    | `VAR_chr1_98039419_C_T;VAR_chr1_98045449_G_C` | If default haplotype: `;` separated list of strings describing events to ignore for haplotype calling in this gene. If non-default haplotype: `;` separated list of strings describing a combination of events that indicates this haplotype. |

For every gene exactly one line should be marked as the default haplotype. This is the haplotype that the reference genome has.

One combination of gene and haplotype name should be marked as the wild-type haplotype for this gene. 
PEACH prefers calling wild-type haplotypes over non-wild-type haplotypes when multiple haplotype combinations are possible.

### Function TSV
| Column    | Example            | Description                           |
|-----------|--------------------|---------------------------------------|
| gene      | `DPYD`             | Name of gene.                         |
| haplotype | `*B3`              | Name of haplotype.                    |
| function  | `Reduced Function` | Description of function of haplotype. |

### Drug TSV
| Column              | Example                                                    | Description                   |
|---------------------|------------------------------------------------------------|-------------------------------|
| gene                | `DPYD`                                                     | Name of gene.                 |
| drugName            | `5-Fluorouracil`                                           | Name of drug.                 |
| urlPrescriptionInfo | `https://www.pharmgkb.org/guidelineAnnotation/PA166104939` | Url with prescription advice. |


## Output
### Events TSV
Name: `[sample_name].peach.events.tsv`

| Column | Example                 | Description                                                        |
|--------|-------------------------|--------------------------------------------------------------------|
| event  | `VAR_chr1_98039419_C_T` | Event ID.                                                          |
| count  | `2`                     | Number of times event was observed. `UNKNOWN` if count is unknown. |

### Gene Events TSV
Name: `[sample_name].peach.gene.events.tsv`

| Column | Example                 | Description                                                        |
|--------|-------------------------|--------------------------------------------------------------------|
| gene   | `DPYD`                  | Name of gene.                                                      |
| event  | `VAR_chr1_98039419_C_T` | Event ID of event linked to gene.                                  |
| count  | `2`                     | Number of times event was observed. `UNKNOWN` if count is unknown. |

### All Haplotypes TSV
Name: `[sample_name].peach.haplotypes.all.tsv`

| Column           | Example  | Description                                                                                                                                                                 |
|------------------|----------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| gene             | `DPYD`   | Name of gene.                                                                                                                                                               |
| combination      | `(*1,2)` | Combination of haplotypes that could explain the observed events for this gene. `;` separated list of pairs of haplotype names and number of times the haplotype is called. |
| count            | `2`      | Number of haplotypes in the combination .                                                                                                                                   |
| nonWildTypeCount | `0`      | Number of non-wild-type haplotypes in the combination .                                                                                                                     |

### Best Haplotypes TSV
Name: `[sample_name].peach.haplotypes.best.tsv`

| Column           | Example                                                          | Description                                                                                                            |
|------------------|------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------|
| gene             | `DPYD`                                                           | Name of gene.                                                                                                          |
| haplotype        | `*1`                                                             | Name of haplotype. `Unresolved Haplotype` if no best haplotype could be determined.                                    |
| count            | `2`                                                              | Number of times this haplotype has been called for this gene.                                                          |
| function         | `Normal Function`                                                | Function of this haplotype. Empty if no `--function_file` argument was provided.                                       |
| linkedDrugs      | `5-Fluorouracil`                                                 | `;` separated list of drugs for which this gene is relevant. Empty if no `--drugs_file` argument was provided.         |
| prescriptionUrls | `https://www.pharmgkb.org/guidelineAnnotation/PA166104939;https` | `;` separated list of urls to prescription advice for related drugs. Empty if no `--drugs_file` argument was provided. |

### QC TSV
Name: `[sample_name].peach.qc.tsv`

| Column           | Example | Description                                                                                                                                                                             |
|------------------|---------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| gene             | `DPYD`  | Name of gene.                                                                                                                                                                           |
| status           | `PASS`  | QC status of calls for this gene. Can be `PASS`, `FAIL_NO_COMBINATION_FOUND`, `FAIL_NO_UNIQUE_BEST_COMBINATION_FOUND`, `FAIL_EVENT_WITH_UNKNOWN_COUNT` or `WARN_TOO_MANY_ALLELES_FOUND` |

## Algorithm
In broad strokes, PEACH does the following:
* Extract relevant calls from the input VCF, where relevance is determined by the configured haplotypes.
* For each gene:
  + Determine for each relevant variant whether the variant is homozygous ref, heterozygous or homozygous alt.
  + Determine the unique simplest combination of haplotypes that completely explains that combination of alt alleles and counts.
    If there is no unique simplest combination of haplotypes that completely explains the combination of alt alleles and counts, then declare "Unresolved Haplotype".
* Create output files.

If no variant call is present in the input VCF for a configured variant, homozygous ref status is assumed.

### Interpret VCF Variant Calls
The calls for the sample `sample_r_id` are extracted from the input VCF, and they are compared to the variants in the config file.
Calls are included when the following are all true:
* Variant is PASS filtered.
* At least one of the covered positions of the call matches at least one of the covered positions of the configured variants.
  In this comparison, the *covered positions* of a call or variant are the positions of the bases in the reference allele.

For each selected call translate the genotype status to a count of occurrences as follows:

| Genotype       | Count   |
|----------------|---------|
| Homozygous ref | 0       |
| Heterozygous   | 1       |
| Homozygous alt | 2       |
| No call        | Unknown |

Note that this assumes that the germline genome is diploid at the relevant sites.

### Determine Valid Haplotype Combinations
Haplotype calling is done for one gene at a time.

First, select all variants for which the following are both true:
* the covered positions overlap with a variant relevant to one of the haplotypes configured for this gene. 
* the variant is not in the list of ignored variants configured for this gene.

If for one of these variant the count is unknown, no haplotypes will be called for this gene.

If all relevant variants have a known count, use recursive descent to determine all combinations of non-default haplotypes that perfectly explain all of these variants and their observed counts.
Non-default, since the default haplotype is the haplotype in the reference genome, and therefore cannot explain any observed variants.
Then, for each haplotype combination with fewer than two called haplotypes, pad the combination with calls of the default haplotype until the combination has two haplotypes.

### Determine the Best Haplotype Combination
Depending on the configured haplotypes, it can be possible for there to be more than one haplotype combination that can explain the observed variants.
As an example, consider the DPYD gene and two variants for that gene: c.1905+1G>A and c.1627A>G.
Separately, these variants form the haplotypes *2A and *5, and the haplotype *2B consists of both of these variants together.
A haplotype can contain multiple variants if these variants have a tendency to be inherited together. 
In this situation both `(*2A,1);(*5,1)` and `(*2B,1);(*1,1)` are valid haplotype combinations.

In the face of such ambiguity PEACH tries to call the "best" haplotype combination. 
First, the preference is for the total number of called haplotypes to be as close as possible to the expected value (2).
Second, wild-type haplotype calls are preferred over non-wild-type haplotypes calls.

The reasoning behind this second preference is that:
* the wild-type haplotype is the wild-type because it is the most common.
* if a haplotype involving multiple variants is configured, this is done because those variants tend to occur together.

In the example above the best haplotype combination would therefore be `(*2B,1);(*1,1)`.

It is not always possible to select a haplotype combination as best. 
When this occurs, the "haplotype" `Unresolved Haplotype` is called instead.

### Peach QC status
PEACH also outputs a QC status per gene. They have the following meaning:

| QC status                               | Meaning                                                                                                         |
|-----------------------------------------|-----------------------------------------------------------------------------------------------------------------|
| `PASS`                                  | No problems were encountered.                                                                                   |
| `WARN_TOO_MANY_ALLELES_FOUND`           | Valid haplotype combinations exist, but all include more than 2 haplotypes calls.                               |
| `FAIL_NO_COMBINATION_FOUND`             | No haplotype combinations could explain all observed variants.                                                  |
| `FAIL_NO_UNIQUE_BEST_COMBINATION_FOUND` | There are multiple valid haplotype combinations, but none of them could be selected as the best.                |
| `FAIL_EVENT_WITH_UNKNOWN_COUNT`         | No haplotype combinations could be called since one of the variants relevant for the gene had an unknown count. |

Any genes with a`FAIL` status get the `Unresolved Haplotype` "haplotype" as the best called haplotype combination.

### Examples
The data in these examples will be the completely fictional.
The examples will focus on fairly "standard" situations, and they will exclude all information that is not necessary to understand these situations.
For details on non-standard situations, see the more detailed subsections of the [Algorithm](#algorithm) section.

Suppose that the config contains the following variants and haplotypes for the fictional gene FAKE,
and that FAKE is the only gene in the config.

| gene | haplotype | default | wildType | events                              |
|------|-----------|---------|----------|-------------------------------------|
| FAKE | *2        | true    | false    | VAR_chr1_400_T_C                    |
| FAKE | *1        | false   | true     | VAR_chr1_100_A_T                    |
| FAKE | *3        | false   | false    | VAR_chr1_100_A_T;VAR_chr1_200_GC_TA |
| FAKE | *4        | false   | false    | VAR_chr1_300_GG_G                   |
| FAKE | *5        | false   | false    | VAR_chr1_100_A_T;VAR_chr1_300_GG_G  |

#### No Calls
If there are no relevant events called, then the only valid haplotype combination that explains these variants is `(*2,2)`.
This is automatically also the best haplotype combination that is called for FAKE.

#### Homozygous Wild-type
Suppose that the called events are the following:

| events             | count |
|--------------------|-------|
| VAR_chr1_100_A_T   | 2     |
| VAR_chr1_200_GC_TA | 0     |
| VAR_chr1_300_GG_G  | 0     |
| VAR_chr1_400_T_C   | 1     |

The call at position 400 is ignored. 
The only haplotype combination that can explain the two variants at position 100 by themselves is `(*1,2)`.
This is automatically also the best haplotype combination that is called for FAKE.

#### Heterozygous Wild-type
Suppose that the called events are the following:

| events             | count |
|--------------------|-------|
| VAR_chr1_100_A_T   | 2     |
| VAR_chr1_200_GC_TA | 1     |
| VAR_chr1_300_GG_G  | 0     |
| VAR_chr1_400_T_C   | 0     |

The only haplotype combination that can explain the two variants at position 100 by themselves is `(*1,1);(*3,1)`.
This is automatically also the best haplotype combination that is called for FAKE.

#### Multiple Valid Haplotypes 1
Suppose that the called events are the following:

| events             | count |
|--------------------|-------|
| VAR_chr1_100_A_T   | 2     |
| VAR_chr1_200_GC_TA | 0     |
| VAR_chr1_300_GG_G  | 2     |
| VAR_chr1_400_T_C   | 0     |

The possible combinations are `(*5,2)`, `(*5,1);(*1,1);(*4,1)` and `(*1,2);(*4,2)`. 
The combination `(*5,2)` is considered to be the best since the total number of haplotypes called is as expected (2).

#### Multiple Valid Haplotypes 2
Suppose that the called events are the following:

| events             | count |
|--------------------|-------|
| VAR_chr1_100_A_T   | 1     |
| VAR_chr1_200_GC_TA | 0     |
| VAR_chr1_300_GG_G  | 1     |
| VAR_chr1_400_T_C   | 0     |

The possible combinations are `(*5,1);(*2;1)` and `(*1,1);(*4,1)`.
In both combinations there are two haplotype calls, so the preference comes down to the number of wild-type haplotypes called.
The best haplotype combination is therefore `(*1,1);(*4,1)`.

## Known Issues / Points for Improvement
* Missing calls are implicitly assumed to represent alleles matching the reference genome. 
This matches the expectation for VCFs from variant calling.
It would be better to remove this implicit assumption and require input VCFs to be genotyped (i.e. contain calls for all relevant positions even when ref).
* Phasing information in the input VCF is currently unused.
* Haplotypes can contain events other than small variants, such as fusions, amplifications and deletions.
These are currently not implemented.
* Germline copy numbers are currently assumed to be 2 for all genes. This is not always true in practice.
* Haplotypes can only be called when all the required events are present the expected number of times. 
Fuzzier matching could potentially better handle uncertainty in the calling in the input files.
* Selection of the best haplotype combination when there are multiple options makes assumptions that don't necessarily hold.
A more careful approach could be beneficial.

## Version History / Download Links
* [2.0.0](https://github.com/hartwigmedical/hmftools/releases/tag/peach-v2.0.0)
  * Converted from Python to Java.
  * Significantly change formats of input and output files.
    * Remove requirement to include both V37 and V38 coordinates in the resource files when calling on a V37 reference genome.
  * Handle overlapping genes correctly.
  * Allow overlapping variants with different ref bases (e.g. deletions and insertions at the same location, overlapping SNVs and MNVs, etc.).
* [1.8](https://github.com/hartwigmedical/peach/releases/tag/v1.8)
  * Update supported Python version from 3.6 to 3.11.
* [1.7](https://github.com/hartwigmedical/peach/releases/tag/v1.7)
  * Adjust PEACH to support UGT1A1 calling.
  * Specifics:
    * Ignore gene of variant in VCF and only use gene names from panel JSON.
    * Allow multiple calls at the same location if reference allele matches.
      * Needed for UGT1A1 *28/*37.
    * Add option to panel JSON to specify variants to ignore.
      * Needed for UGT1A1 *36.
    * Stop treating empty variant annotation as unknown variant annotation, but instead indicate as annotation="NONE".
* [1.6](https://github.com/hartwigmedical/peach/releases/tag/v1.6)
  * Make --sample_t_id argument optional.
  * Fix crash when PAVE_TI has "Number=." in VCF header.
* [1.5](https://github.com/hartwigmedical/peach/releases/tag/v1.5)
  * Add support for [PAVE](https://github.com/hartwigmedical/hmftools/tree/master/pave) annotations.
  * Add optional "canonicalTranscript" entry in panel JSON.
    * This entry is required for parsing PAVE annotations, and ignored for SnpEff annotations.
* [1.4](https://github.com/hartwigmedical/peach/releases/tag/v1.4)
  * Change formats of output files to essentially being PEACH v1.0 output files with some additional columns,
    to avoid breaking the expectations from downstream tools based on semantic versioning.
* [1.3](https://github.com/hartwigmedical/peach/releases/tag/v1.3)
  * Update scikit-allele version to fix pip-install failure.
* [1.2](https://github.com/hartwigmedical/peach/releases/tag/v1.2)
  * Allow for different chromosome names wrt v37 and v38
  * Includes changes to the panel JSON and the genotype output file.
  * Improve clarity and consistency of logging.
* [1.1](https://github.com/hartwigmedical/peach/releases/tag/v1.1)
  * Add shell script `peach` for running PEACH.
  * Remove VCF filtering step.
    * Remove VCFTools dependency.
    * Remove filtered-VCF output file.
  * Add experimental support for input VCF's for reference genomes with version v38.
  * Change format of arguments to PEACH.
    * Arguments are no longer positional.
    * Remove arguments vcftools, recreate_bed and transcript_tsv, since they are no longer needed.
    * Add optional parameter for experimental v38 VCF support.
  * Change format of panel JSON.
    * Change key "url_prescription_info" to "urlPrescriptionInfo" for consistency with other keys.
    * Add "annotationV37" key for reference sequence differences, for support of v38 reference genomes.
  * Adjust format of genotype TSV output file.
    * Split "haplotype" column into "haplotype" and "zygosity".
  * Add script for running tests.
* [1.0](https://github.com/hartwigmedical/peach/releases/tag/v1.0)
  * First release
