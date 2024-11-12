CHORD: Classifier of HOmologous Recombination Deficiency
================

CHORD is a random forest model that predicts homologous recombination deficiency (HRD) using somatic mutations. CHORD primarily uses the 
relative counts of somatic mutation contexts, namely deletions with flanking microhomology and 1-100kb structural duplications, with 
1-100kb structural duplications being used to distinguish BRCA1-type HRD from BRCA2-type HRD.

This repo contains Java the implementation of the [mutSigExtractor](https://github.com/UMCUGenetics/mutSigExtractor) and 
[CHORD](https://github.com/UMCUGenetics/CHORD) R packages.

For more info on CHORD, please see https://www.nature.com/articles/s41467-020-19406-4: </br> 
**Pan-cancer landscape of homologous recombination deficiency**. </br>
*Luan Nguyen, John Martens, Arne Van Hoeck, Edwin Cuppen. Nat Commun 11, 5584 (2020).*

# Usage

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
`-snv_indel_vcf_file` and `-sv_vcf_file`, 

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

We can run the **feature extraction** step by calling `ChordDataPrep`. This Java class replaces the functionality of the 
[mutSigExtractor](https://github.com/UMCUGenetics/mutSigExtractor) R package.

```shell
java -cp chord.jar com.hartwig.hmftools.chord.prep.ChordDataPrep \
    -sample SAMPLE_1,SAMPLE_2 \
    -snv_indel_vcf_file "/data/datasets/*/purple/*.purple.somatic.vcf.gz" \
    -sv_vcf_file "/data/datasets/*/purple/*.purple.sv.vcf.gz" \
    -output_dir /path/to/output/dir/ \
    -ref_genome Homo_sapiens.GRCh37.GATK.illumina.fasta \
```

We can run the random forest **prediction** step by calling `ChordModel`. This Java class replaces the 
[CHORD](https://github.com/UMCUGenetics/CHORD) R package, and is a wrapper around a simple R script to run the CHORD random forest.
```shell
java -cp chord.jar com.hartwig.hmftools.chord.predict.ChordModel \
    -mut_contexts_file TUMOR_SAMPLE.chord.mutation_contexts.modified.tsv \
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
         A[C>A]A A[C>A]C ... del.mh.bimh.1 del.mh.bimh.2 ... DEL_0e00_1e03_bp DEL_1e03_1e04_bp ...
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
