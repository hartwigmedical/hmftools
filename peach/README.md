# PEACH

**P**harmacogenomic **E**valuator **A**nd **C**aller of **H**aplotypes (PEACH) is a pharmacogenomics tool developed for
the [Hartwig Medical Foundation pipeline](https://github.com/hartwigmedical/pipeline5).
It imports haplotypes and related variants from config files, reports the presence of these variants in a
germline VCF, and infers the simplest combination of haplotypes that explains the presence of these variants.

An earlier Python version of this tool is available [here](https://github.com/hartwigmedical/peach).

## Contents
TODO: fix this later

* [Installation](#installation)
* [Arguments](#arguments)
  + [Mandatory Arguments](#mandatory-arguments)
  + [Optional Arguments](#optional-arguments)
* [Input](#input)
  + [VCF](#vcf)
  + [Config files](#config-files)
* [Output](#output)
  + [Genotype TSV file](#genotype-tsv-file)
  + [Calls TSV file](#calls-tsv-file)
* [Algorithm](#algorithm)
  + [Get_VCF Variant Calls](#get-vcf-variant-calls)
  + [Infer Haplotypes](#infer-haplotypes)
  + [Examples](#examples)
* [Known issues / points for improvement](#known-issues--points-for-improvement)
* [Version History / Download Links](#version-history--download-links)

## Installation
To install, download the latest compiled jar file from the [download links](#version-history-and-download-links).

PEACH requires Java 11+.

TODO: Config files

## Arguments
#### Example Usage
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
| wildType  | `false`                                       | `true` if this haplotype is considered to be the wild type haplotype of this gene. `false` otherwise.                                                                                                                                         |
| events    | `VAR_chr1_98039419_C_T;VAR_chr1_98045449_G_C` | If default haplotype: `;` separated list of strings describing events to ignore for haplotype calling in this gene. If non-default haplotype: `;` separated list of strings describing a combination of events that indicates this haplotype. |

For every gene exactly one line should be marked as the default haplotype. This is the haplotype that the reference genome has.

One combination of gene and haplotype name should be marked as the wild type haplotype for this gene.

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

| Column | Example                 | Description                         |
|--------|-------------------------|-------------------------------------|
| event  | `VAR_chr1_98039419_C_T` | Event ID.                           |
| count  | `2`                     | Number of times event was observed. |

### Gene Events TSV
Name: `[sample_name].peach.gene.events.tsv`

| Column | Example                 | Description                         |
|--------|-------------------------|-------------------------------------|
| gene   | `DPYD`                  | Name of gene.                       |
| event  | `VAR_chr1_98039419_C_T` | Event ID of event linked to gene.   |
| count  | `2`                     | Number of times event was observed. |

### All Haplotypes TSV
Name: `[sample_name].peach.haplotypes.all.tsv`

| Column           | Example  | Description                                                                                                                                                                 |
|------------------|----------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| gene             | `DPYD`   | Name of gene.                                                                                                                                                               |
| combination      | `(*1,2)` | Combination of haplotypes that could explain the observed events for this gene. `;` separated list of pairs of haplotype names and number of times the haplotype is called. |
| nonWildTypeCount | `0`      | Number of non wild type haplotypes in the combination .                                                                                                                     |

### Best Haplotypes TSV
Name: `[sample_name].peach.haplotypes.best.tsv`

| Column           | Example                                                          | Description                                                                                                            |
|------------------|------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------|
| gene             | `DPYD`                                                           | Name of gene.                                                                                                          |
| haplotype        | `*1`                                                             | Name of haplotype.                                                                                                     |
| count            | `2`                                                              | Number of times this haplotype has been called for this gene.                                                          |
| function         | `Normal Function`                                                | Function of this haplotype. Empty if no `--function_file` argument was provided.                                       |
| linkedDrugs      | `5-Fluorouracil`                                                 | `;` separated list of drugs for which this gene is relevant. Empty if no `--drugs_file` argument was provided.         |
| prescriptionUrls | `https://www.pharmgkb.org/guidelineAnnotation/PA166104939;https` | `;` separated list of urls to prescription advice for related drugs. Empty if no `--drugs_file` argument was provided. |

### QC TSV
Name: `[sample_name].peach.qc.tsv`

| Column           | Example | Description                                                                                                                                            |
|------------------|---------|--------------------------------------------------------------------------------------------------------------------------------------------------------|
| gene             | `DPYD`  | Name of gene.                                                                                                                                          |
| status           | `PASS`  | QC status of calls for this gene. Can be `PASS`, `FAIL_NO_COMBINATION_FOUND`, `FAIL_NO_UNIQUE_BEST_COMBINATION_FOUND` or `WARN_TOO_MANY_ALLELES_FOUND` |

## Algorithm
TODO: how to handle 37 calls needs to go somewhere

TODO: describe ignored variants

In broad strokes, PEACH does the following:
* Extract relevant calls from VCF, where relevance is determined by the configured haplotypes.
* For each gene:
  + Determine for each relevant variant whether the variant is homozygous ref, heterozygous or homozygous alt.
  + Determine the unique simplest combination of haplotypes that completely explains that combination of alt alleles and counts.
    If there is no unique simplest combination of haplotypes that completely explains the combination of alt alleles and counts, then declare "Unresolved Haplotype".
* Create output files.

If no variant call is present in the input VCF for a configured variant, homozygous ref is assumed.

### Get VCF Variant Calls
TODO: clean this up

The calls for the sample `sample_r_id` are extracted from the input VCF, and they are compared to the variants in the config file.
Calls are included when the following are all true:
* Variant is PASS filtered.
* At least one of the covered positions of the call matches at least one of the covered positions of the configured variants.
  In this comparison, the *covered positions* of a call or variant are the positions of the bases in the reference allele.



### Infer Haplotypes

TODO: be consistent and clear about the difference between default and wild-type allele.

TODO: Properly explain what wild type is

The goal is to find the simplest combination of haplotypes that explains the called variants.

Sometimes, more than one combination of haplotypes could explain the calls.
As an example, consider the DPYD gene and two variants for that gene: c.1905+1G>A and c.1627A>G.
Separately, these variants form the haplotypes *2A and *5, and the haplotype *2B consists of both of these variants together.
A haplotype can contain multiple variants if these variants have a tendency to be inherited together.
If each variant is called once and if all of these variants and haplotypes are included in the variant configuration, then, 
to take this combined inheritance into account,
PEACH prefers to call the haplotype combination as *2B_HET/*1_HET and not as *2A_HET/*5_HET.
If you want PEACH to call *2A_HET/*5_HET instead of *2B_HET/*1_HET in this situation, then simply don't include *2B in the configuration.

To make this more precise, define the *length* of a haplotype combination as the total number of non-wild-type
haplotypes in the combination, where homozygous haplotype calls are counted as 2.
PEACH will always attempt to call the unique haplotype combination of minimum length that explain all the variant calls.
If there are no haplotype combinations that explain all the variant calls,
or if there is more than one combination of the same minimum length,
then the haplotype combination for that gene is called as "Unresolved Haplotype".

The only valid haplotype combination of length 0 is the homozygous default haplotype.
Valid haplotype combinations of length 1 always include precisely one heterozygous default haplotype call.
Valid haplotype combinations of length at least 2 do not contain any calls for the default haplotype.

Note that when at least one of the VCF calls overlaps with but is not identical to
one of the variants in the configuration, then the haplotype combination "Unresolved Haplotype" will be called,
because this variant will be an event that is not part of any haplotypes in the configuration.

#### Haplotype Calling Algorithm

TODO: introduce and consistently use the concept of "best haplotype combination"

TODO: Mention assumption of 2 haplotypes per gene.

Haplotypes are called for each gene separately. First, collect the calls that are included in at least one non-default haplotype.
Use recursive descent to determine all combinations of non-default haplotypes that perfectly explain all of these variants.
If there are no such combinations, then no haplotype combination can be called for this gene.
If such combinations do exist, then the next step is to determine the length of
each valid haplotype combination. Find the minimum length of the valid haplotype combinations,
and select the haplotype combinations whose length is equal to this minimum.
If precisely one such haplotype combination exists, then this combination will be called for this gene.
If more than one haplotype combination of minimum length exists, then no haplotype combination is called for this gene.

#### Peach QC status

TODO

### Examples
TODO: Fix example

The data in these examples will be the completely fictional.
The examples will focus on fairly "standard" situations, and they will exclude all information that is not necessary to understand these situations.
For details on non-standard situations, see the more detailed subsections of the [Algorithm](#algorithm) section.

Suppose that the config contains the following variants and haplotypes for the fictional gene FAKE,
and that FAKE is the only gene in the config.

| Rs Id | Reference Allele V37 | Reference Allele V38 | V38 Annotation for Reference Sequence Difference |
|-------|----------------------|----------------------|--------------------------------------------------|
| rs1   | A                    | A                    | N/A                                              |
| rs2   | TA                   | GC                   | c.6543GC>TA                                      |
| rs3   | GG                   | GG                   | N/A                                              |

| Haplotype      | Variants (Rs Id: Variant Allele wrt v38) |
|----------------|------------------------------------------|
| *1 (wild type) | None                                     |
| *2             | rs1: T                                   |
| *3             | rs2: TA                                  |
| *4             | rs3: G                                   |
| *5             | rs1: T, rs3: G                           |

#### No Calls
If there are no calls wrt v37 in the VCF, then the dual calls are:

| Rs Id | Allele1 | Allele2 | Variant Annotation V37 | Filter V37 | Variant Annotation V38 | Filter V38    |
|-------|---------|---------|------------------------|------------|------------------------|---------------|
| rs1   | A       | A       | REF_CALL               | NO_CALL    | REF_CALL               | NO_CALL       |
| rs2   | TA      | TA      | REF_CALL               | NO_CALL    | c.6543GC>TA            | INFERRED_PASS |
| rs3   | GG      | GG      | REF_CALL               | NO_CALL    | REF_CALL               | NO_CALL       |

The only valid haplotype combination that explains these variants is *3_HOM,
so this is the haplotype combination that is called for FAKE.

#### Homozygous Wild Type
Suppose that the v37 calls are the following:

| Rs Id | Allele1 | Allele2 | Variant Annotation V37 | Filter V37 |
|-------|---------|---------|------------------------|------------|
| rs1   | A       | A       | c.8483A>T              | PASS       |
| rs2   | GC      | GC      | c.6543TA>GC            | PASS       |
| rs3   | GG      | GG      | c.4838GG>G             | PASS       |

In this case, the dual calls are:

| Rs Id | Allele1 | Allele2 | Variant Annotation V37 | Filter V37 | Variant Annotation V38 | Filter V38 |
|-------|---------|---------|------------------------|------------|------------------------|------------|
| rs1   | A       | A       | REF_CALL               | PASS       | REF_CALL               | PASS       |
| rs2   | GC      | GC      | c.6543TA>GC            | PASS       | REF_CALL               | PASS       |
| rs3   | GG      | GG      | REF_CALL               | PASS       | REF_CALL               | PASS       |

The only valid haplotype combination is *1_HOM, so this haplotype combination is called for FAKE.

#### Heterozygous Wild Type
Suppose that the v37 calls are the following:

| Rs Id | Allele1 | Allele2 | Variant Annotation V37 | Filter V37 |
|-------|---------|---------|------------------------|------------|
| rs1   | A       | A       | c.8483A>T              | PASS       |
| rs2   | TA      | GC      | c.6543TA>GC            | PASS       |
| rs3   | GG      | GG      | c.4838GG>G             | PASS       |

The dual calls are:

| Rs Id | Allele1 | Allele2 | Variant Annotation V37 | Filter V37 | Variant Annotation V38 | Filter V38 |
|-------|---------|---------|------------------------|------------|------------------------|------------|
| rs1   | A       | A       | REF_CALL               | PASS       | REF_CALL               | PASS       |
| rs2   | GC      | TA      | c.6543TA>GC            | PASS       | c.6543GC>TA            | PASS       |
| rs3   | GG      | GG      | REF_CALL               | PASS       | REF_CALL               | PASS       |

The only valid haplotype combination for FAKE is *3_HET/*1_HET, so this haplotype combination is called.

#### Multiple Valid Haplotypes
Suppose that the v37 calls are the following:

| Rs Id | Allele1 | Allele2 | Variant Annotation V37 | Filter V37 |
|-------|---------|---------|------------------------|------------|
| rs1   | A       | T       | c.8483A>T              | PASS       |
| rs2   | GC      | GC      | c.6543TA>GC            | PASS       |
| rs3   | G       | G       | c.4838GG>G             | PASS       |

The resulting dual calls are:

| Rs Id | Allele1 | Allele2 | Variant Annotation V37 | Filter V37 | Variant Annotation V38 | Filter V38 |
|-------|---------|---------|------------------------|------------|------------------------|------------|
| rs1   | A       | T       | c.8483A>T              | PASS       | c.8483A>T              | PASS       |
| rs2   | GC      | GC      | c.6543TA>GC            | PASS       | REF_CALL               | PASS       |
| rs3   | G       | G       | c.4838GG>G             | PASS       | c.4838GG>G             | PASS       |

The valid haplotype combinations are *2_HET/*4_HOM and *4_HET/*5_HET.
These combinations have lengths 3 and 2, respectively, so the second combination is preferred.
The called haplotype combination for FAKE is *4_HET/*5_HET.

## Known issues / points for improvement
TODO: Mention genotype vs variant calling
TODO: Use phasing
TODO: Other event types
TODO: Take different germline copy numbers into consideration
TODO: Try to work on tumor-only
TODO: Improve handling of uncertainty and unknown haplotypes

## Version History / Download Links
* Upcoming
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
