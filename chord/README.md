CHORD: Classifier of HOmologous Recombination Deficiency
================

CHORD is a random forest model that predicts homologous recombination deficiency (HRD) using relative counts of somatic mutation contexts, 
primarily deletions with flanking microhomology to distinguish HRD vs non-HRD, and 1-100kb duplications to distinguish BRCA1-type vs 
BRCA2-type HRD. For more info on CHORD, please see the paper: [Pan-cancer landscape of homologous recombination deficiency](https://www.nature.com/articles/s41467-020-19406-4).

`hmftools/chord` is a Java reimplementation of the original two R packages used to run CHORD:
- [mutSigExtractor](https://github.com/UMCUGenetics/mutSigExtractor): Performs the feature extraction from VCFs. The core functionality has been entirely migrated to Java
- [CHORD](https://github.com/UMCUGenetics/CHORD): Runs the CHORD random forest. This is now a simple R script with a Java wrapper

# Usage

## Arguments

**Main class** (`com.hartwig.hmftools.chord.ChordApplication`) **and feature extraction** (`com.hartwig.hmftools.chord.prep.ChordDataPrep`)

| Argument                | Example                                           | Description                                                                                                       |
|-------------------------|---------------------------------------------------|-------------------------------------------------------------------------------------------------------------------|
| `-sample`               | TUMOR_SAMPLE                                      | <sup>1</sup> Sample name                                                                                          |
| `-sample_id_file`       | sample_ids.txt                                    | <sup>1</sup> A one-column text file listing sample IDs                                                            |
| `-snv_indel_vcf_file`   | "/data/datasets/*/purple/*.purple.somatic.vcf.gz" | <sup>2</sup> Path to a VCF containing SNVs and INDELs                                                             |
| `-sv_vcf_file`          | "/data/datasets/*/purple/*.purple.sv.vcf.gz"      | <sup>2</sup> Path to a VCF containing SVs                                                                         |
| `-purple_dir`           | "/data/datasets/*/purple/"                        | <sup>2</sup> Directory containing the PURPLE files                                                                |
| `-output_dir`           | output/                                           | Directory to write the output files                                                                               |
| `-ref_genome`           | Homo_sapiens.GRCh37.GATK.illumina.fasta           | Path to a reference genome .fasta file. The .dict and .fai index files must also be present in the same directory |
| `-include_non_pass`     |                                                   | **Flag**. Include non PASS variants during feature extraction                                                     |
| `-write_detailed_files` |                                                   | **Flag**. Write TSV files containing per-sample and per-variant feature extraction information                    |
| `-threads`              | 8                                                 | Number of threads to use. Each thread processes one sample at a time                                              |
| `-log_level`            | DEBUG                                             | Set log level to one of: ERROR, WARN, INFO, DEBUG or TRACE                                                        |
| `-log_debug`            |                                                   | **Flag**. Set log level to DEBUG                                                                                  |

Notes:
1. Provide either `-sample` or `-sample_id_file`
2. Provide either `-snv_indel_vcf_file` and `-sv_vcf_file` together, or `-purple_dir`

**Predicting from existing mutation contexts** (`com.hartwig.hmftools.chord.predict.ChordModel`)

| Argument             | Example                                  | Description                                                   |
|----------------------|------------------------------------------|---------------------------------------------------------------|
| `-mut_contexts_file` | TUMOR_SAMPLE.chord.mutation_contexts.tsv | Path to the mutation contexts file created by `ChordDataPrep` |
| `-output_file`       | TUMOR_SAMPLE.chord.prediction.tsv        | Path output the predictions                                   |


## Single sample mode

CHORD can be run using the below example command:

```shell
java -jar chord.jar \
    -sample TUMOR_SAMPLE \
    -snv_indel_vcf_file TUMOR_SAMPLE.purple.somatic.vcf.gz \
    -sv_vcf_file TUMOR_SAMPLE.purple.sv.vcf.gz \
    -output_dir /path/to/output/dir/ \
    -ref_genome Homo_sapiens.GRCh37.GATK.illumina.fasta
```

> [!NOTE]
> Make sure that the provided reference genome matches the version that the VCF files were generated with.

## Multi-sample mode

To run CHORD in multi-sample mode, we can provide comma-delimited sample IDs to `-sample`, and paths with wildcards (`*`) to 
`-snv_indel_vcf_file` and `-sv_vcf_file`. We can also provide `-threads` to run in multithreaded mode.

```shell
java -jar chord.jar \
    -sample SAMPLE_1,SAMPLE_2 \
    -snv_indel_vcf_file "/data/datasets/*/purple/*.purple.somatic.vcf.gz" \
    -sv_vcf_file "/data/datasets/*/purple/*.purple.sv.vcf.gz" \
    -output_dir /path/to/output/dir/ \
    -ref_genome Homo_sapiens.GRCh37.GATK.illumina.fasta \
    -threads 2
```

> [!NOTE]
> Paths with wildcards must be surrounded by double quotes (`"`)

The above command internally generates the following paths:
```
/data/datasets/SAMPLE_1/purple/SAMPLE_1.purple.somatic.vcf.gz
/data/datasets/SAMPLE_1/purple/SAMPLE_1.purple.sv.vcf.gz
/data/datasets/SAMPLE_2/purple/SAMPLE_2.purple.somatic.vcf.gz
/data/datasets/SAMPLE_2/purple/SAMPLE_2.purple.sv.vcf.gz
```

## Alternative input arguments
As shown in the below command, we can provide:
- `-sample_id_file` instead of listing sample IDs to `sample`
- `-purple_dir` instead of `-snv_indel_vcf_file` and `-sv_vcf_file`

```shell
java -jar chord.jar \
    -sample_id_file sample_ids.txt \
    -purple_dir "/data/datasets/*/purple" \
    -output_dir /path/to/output/dir/ \
    -ref_genome Homo_sapiens.GRCh37.GATK.illumina.fasta \
    -threads 2
```

The `-sample_id_file` is a one-column text file listing the sample IDs:
```
SampleId
SAMPLE_1
SAMPLE_2
```

## Running feature extraction or prediction separately

We can run the **feature extraction** step by calling `ChordDataPrep`:

```shell
java -cp chord.jar com.hartwig.hmftools.chord.prep.ChordDataPrep \
    -sample SAMPLE_1,SAMPLE_2 \
    -snv_indel_vcf_file "/data/datasets/*/purple/*.purple.somatic.vcf.gz" \
    -sv_vcf_file "/data/datasets/*/purple/*.purple.sv.vcf.gz" \
    -output_dir /path/to/output/dir/ \
    -ref_genome Homo_sapiens.GRCh37.GATK.illumina.fasta \
```

We can run the random forest **prediction** step by calling `ChordModel`:
```shell
java -cp chord.jar com.hartwig.hmftools.chord.predict.ChordModel \
    -mut_contexts_file TUMOR_SAMPLE.chord.mutation_contexts.tsv \
    -output_file TUMOR_SAMPLE.chord.prediction.tsv 
```

# Output

CHORD will output 2 files:
- `*.chord.mutation_contexts.tsv`: Input features to CHORD (SNV, indel and SV context counts)
- `*.chord.prediction.tsv`: Predictions from the CHORD random forest

These files can contain the data for one or more samples.

## Mutation contexts

The `*.chord.mutation_contexts.tsv` file contains the counts of SNV, INDEL and SV contexts:

```
sample_id A[C>A]A A[C>A]C ... del.mh.bimh.1 del.mh.bimh.2 ... DEL_0e00_1e03_bp DEL_1e03_1e04_bp ...
 SAMPLE_1      76      41 ...             5            16 ...                5                8 ...
 SAMPLE_2      66      63 ...            19            16 ...                3                7 ...
```

The mutation contexts include:
- SNV trinucleotide contexts, i.e. split by:
  - Substitution type
  - 5' and 3' flanking base
- Insertions and deletions (`del`,`ins`) split by type:
  - Repeats split by repeat length (`rep.len`)
  - Microhomology split by number of bases in microhomology (`mh.bimh`)
  - Other indels split by sequence length (`none.len`)
- SVs split by:
  - Type: `DEL`,`DUP`,`INV`,`TRA`
  - Length bins with intervals: 0, 1e3, 1e4, 1e5, 1e6, 1e7, infinity

## Predictions

The below example `*.chord.prediction.tsv` file shows different possible predictions that CHORD can produce: 

```
  sample  p_BRCA1  p_BRCA2  p_hrd             hr_status              hrd_type  remarks_hr_status  remarks_hrd_type
SAMPLE_1    0.000    0.000  0.000         HR_proficient                  none  
SAMPLE_2    0.798    0.056  0.854          HR_deficient            BRCA1_type  
SAMPLE_3    0.062    0.680  0.742          HR_deficient            BRCA2_type  
SAMPLE_4    0.000    0.002  0.002  cannot_be_determined                  none         <50 indels
SAMPLE_5    0.331    0.410  0.741          HR_deficient  cannot_be_determined                              <30 SVs
```

### Probabilities and HR status
- `p_hrd`: Probability of HRD
- `p_BRCA1`: Probability of BRCA1-type HRD
- `p_BRCA2`: Probability of BRCA2-type HRD
- `hr_status`: Indicates if a sample is `HR_deficient` (`p_hrd` \>= 0.5) or `HR_proficient`.
- `hrd_type`: If `hr_status` is `HR_deficient`, indicates if the sample has `BRCA1-type` or `BRCA2-type` HRD; i.e. max(`p_BRCA1`,`p_BRCA2`)

### Quality control checks
- CHORD requires **\>=100 indels** to accurately determine whether a sample is HRD. If this criterion is not met, `hr_status` will be
`cannot_be_determined` and `remarks_hr_status` will be `<50 indels`.
- CHORD **cannot be applied to MSI samples**. If an MSI sample is detected, `hr_status` will be `cannot_be_determined` and
`remarks_hr_status` will be `Has MSI (>14000 indel.rep)`
- CHORD requires **\>=30 SVs** to accurately determine HRD subtype. If this criterion is not met, `hrd_type` will be
`cannot_be_determined`, and `remarks_hr_status` will be `<30 SVs`.

The user may of course ignore these remarks and proceed with using the raw probabilities outputted by CHORD (`p_hrd` and/or`p_BRCA1`/`pBRCA2`).
