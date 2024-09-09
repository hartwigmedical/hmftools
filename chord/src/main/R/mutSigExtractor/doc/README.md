mutSigExtractor
================

## Description

mutSigExtractor is an R package for extracting mutation contexts from
vcf files. Extraction can be performed for the following mutation types:

  - **SNV**: 96 trinucleotide contexts, consisting of the point mutation
    and the 5’ and 3’ flanking nucleotides
  - **DBS**: 78 possible double base substitution types
  - **Indel**: Includes indels within repeat regions, indels with
    flanking microhomology; and other indels. Each category is further
    stratified by the repeat unit length, the number of bases in the
    indel sequence that are homologous, and the indel sequence length,
    respectively. Alternatively, the PCAWG indel contexts can also be
    extracted
  - **SV**: structural variants stratified by type (deletions,
    duplications, inversions, translocations) and length (0 to \>10Mb)

Signatures can also be extracted from these contexts for SNVs, indels,
DBSs (COSMIC/PCAWG), as well as for SVs ([Nik-Zainal et
al. 2016](https://www.nature.com/articles/nature17676)).

## Installation

mutSigExtractor requires some bioconductor packages to first be
installed.

``` r
## Bioconductor packages required by mutSigExtractor
install.packages('BiocManager')
BiocManager::install('BSgenome') ## Install genome parser
BiocManager::install('BSgenome.Hsapiens.UCSC.hg19') ## Install the default genome
BiocManager::install('GenomeInfoDb')

## Install mutSigExtractor directly from github using devtools
install.packages("devtools")
devtools::install_github('https://github.com/UMCUGenetics/mutSigExtractor/')
```

Other genomes (e.g. for hg38: `BSgenome.Hsapiens.UCSC.hg38`) can also be
used. Please see the below tutorial for details.

## Tutorial

The `COLO829v003T.purple.somatic.vcf.gz` and
`COLO829v003T.purple.sv.vcf.gz` files in `doc/vcf/` will be used to
demonstrate how to use mutSigExtractor.

``` r
setwd('/path/to/mutSigExtractor')
vcf_snv <- 'doc/vcf/COLO829v003T.purple.somatic.vcf.gz'
vcf_sv <- 'doc/vcf/COLO829v003T.purple.sv.vcf.gz'
```

The main functions for extracting signatures are:

    extractSigsSnv()
    extractSigsDbs()
    extractSigsIndel()
    extractSigsSv()

Note that SNVs, DBSs and indels are often reported in the same vcf file.
Therefore, `extractSigsSnv()`,`extractSigsDbs()`, and
`extractSigsIndel()` will automatically select the relevant mutation
types.

### Extracting contexts and signatures

#### SNVs

Here, extraction of SNV contexts and signatures on one sample using
`extractSigsSnv()` will be demonstrated. The same concepts shown here
can be applied to the other mutation type functions.

With the below code we can extract the 96 trinucleotide contexts. This
returns a single column matrix of mutation counts per context.

Why is the output not simply a vector? This is so that the output can be
written to a txt file, which is useful when processing a large number of
samples on an HPC. When processing multiple samples locally, one can
simply use `cbind()` combine the context counts from all samples into
one matrix.

Note that, it is recommended that the `vcf.filter` argument is set to
‘PASS’ (or ‘.’ for certain vcf files) to remove low quality variants.

``` r
contexts_snv <- extractSigsSnv(vcf.file=vcf_snv, vcf.filter='PASS', output='contexts')
head(contexts_snv)
```

    ##         COLO829v003T.purple.somatic.vcf.gz
    ## A[C>A]A                                132
    ## A[C>A]C                                 62
    ## A[C>A]G                                 17
    ## A[C>A]T                                 57
    ## C[C>A]A                               1892
    ## C[C>A]C                                237

To extract signatures, we can then fit these contexts to e.g. the COSMIC
or PCAWG signature profiles. These are included in the package as
`SBS_SIGNATURE_PROFILES_V2` and `SBS_SIGNATURE_PROFILES_V3`
respectively. For this example we will use the PCAWG profiles
(`SBS_SIGNATURE_PROFILES_V3`).

``` r
SBS_SIGNATURE_PROFILES_V3[1:5,1:5]
```

    ##             SBS1     SBS2    SBS3   SBS4    SBS5
    ## A[C>A]A 0.000886 5.80e-07 0.02080 0.0422 0.01200
    ## A[C>A]C 0.002280 1.48e-04 0.01650 0.0333 0.00944
    ## A[C>A]G 0.000177 5.23e-05 0.00175 0.0156 0.00185
    ## A[C>A]T 0.001280 9.78e-05 0.01220 0.0295 0.00661
    ## C[C>A]A 0.000312 2.08e-04 0.02250 0.0807 0.00743

Signature fitting is done using `fitToSignatures()` This function uses
the non-negative linear least squares algorithm based on `lsqnonneq()`
from the `pracma` package. `mut.context.counts` can be a numeric vector,
matrix, or dataframe. If a matrix/dataframe, rows represent samples and
columns represent contexts.

``` r
sigs_snv <- fitToSignatures(
   mut.context.counts=contexts_snv[,1], 
   signature.profiles=SBS_SIGNATURE_PROFILES_V3
)
head(sigs_snv)
```

    ##      SBS1      SBS2      SBS3      SBS4      SBS5      SBS6 
    ##  71.02426 770.66129   0.00000   0.00000   0.00000   0.00000

Alternatively, signatures can be extracted directly from the vcf by
specifying `output='signatures'` in `extractSigsSnv()`

``` r
sig_snv_2 <- extractSigsSnv(
   vcf.file=vcf_snv, vcf.filter='PASS', output='signatures', 
   signature.profiles=SBS_SIGNATURE_PROFILES_V3
)
head(sig_snv_2)
```

    ##      COLO829v003T.purple.somatic.vcf.gz
    ## SBS1                           71.02426
    ## SBS2                          770.66129
    ## SBS3                            0.00000
    ## SBS4                            0.00000
    ## SBS5                            0.00000
    ## SBS6                            0.00000

### Output contexts for DBSs, indels, and SVs

The context extractions for mutation types other than SNVs are shown
below.

#### DBSs

``` r
contexts_dbs <- extractSigsDbs(vcf.file=vcf_snv, vcf.filter='PASS', output='contexts')
head(contexts_dbs)
```

    ##       COLO829v003T.purple.somatic.vcf.gz
    ## AC>CA                                  0
    ## AC>CG                                  0
    ## AC>CT                                  5
    ## AC>GA                                  3
    ## AC>GG                                  0
    ## AC>GT                                  1

#### Indels

For indels, `extractSigsIndel()` defaults to `method='CHORD'` which is
used by Classifier of HOmologous Recombination Deficiency (CHORD). When
`method='CHORD'`, contexts are always extracted (i.e. no `output`
argument).

``` r
contexts_indel <- extractSigsIndel(vcf.file=vcf_snv, vcf.filter='PASS')
head(contexts_indel)
```

    ##               COLO829v003T.purple.somatic.vcf.gz
    ## del.rep.len.1                                155
    ## del.rep.len.2                                 24
    ## del.rep.len.3                                 11
    ## del.rep.len.4                                 13
    ## del.rep.len.5                                 24
    ## ins.rep.len.1                                121

However, the PCAWG indel contexts can also be extracted by setting
`method='PCAWG'`. Here, `output` can be `'signatures'` to directly
extract the PCAWG indel signatures.

``` r
contexts_indel <- extractSigsIndel(vcf.file=vcf_snv, vcf.filter='PASS', method='PCAWG', output='contexts')
head(contexts_indel)
```

    ##            COLO829v003T.purple.somatic.vcf.gz
    ## del.1.C.1                                  25
    ## del.1.C.2                                  20
    ## del.1.C.3                                   8
    ## del.1.C.4                                   6
    ## del.1.C.5                                   1
    ## del.1.C.6+                                  4

#### SVs

SV vcf files generally do not adhere to one standard. `extractSigsSv()`
currently supports SV vcf parsing of GRIDSS (conforms to vcf spec 4.2)
and Manta vcfs. You can specify the vcf type with `sv.caller='gridss'`
(default) or `sv.caller='manta'`.

``` r
contexts_sv <- extractSigsSv(vcf.file=vcf_sv, vcf.filter='PASS', output='contexts', sv.caller='gridss')
head(contexts_sv)
```

    ##                  COLO829v003T.purple.sv.vcf.gz
    ## DEL_0e00_1e03_bp                             6
    ## DEL_1e03_1e04_bp                             8
    ## DEL_1e04_1e05_bp                            11
    ## DEL_1e05_1e06_bp                            12
    ## DEL_1e06_1e07_bp                             1
    ## DEL_1e07_Inf_bp                              0

In case you want to use SV vcfs from other callers, you can provide a
dataframe as input to `extractSigsSv()` containing SV type and length
info. See below for more info.

### Dataframes as input

Alternatively, dataframes can be used as input, which is handy if you
want to parse the vcfs yourself, or have inputs in tabular format.

For SNVs and indels, a dataframe with the columns: chrom, pos, ref, alt.

    ##   CHROM      POS REF ALT
    ## 1     1 16145827   A   C
    ## 2     1 16492085   G   C
    ## 3     1 17890303   C   G
    ## 4     1 18877885   G   A
    ## 5     1 18919776   T   C

For SVs, a dataframe with the columns: sv\_type, sv\_len. For sv\_type,
values must be DEL, DUP, INV, TRA (deletions, duplications, inversions,
translocations). For translocations, sv\_len information is discarded.
The column names themselves do not matter, as long as the columns are in
the aforementioned order.

    ##   sv_type    sv_len
    ## 1     TRA        NA
    ## 2     DEL      1696
    ## 3     DEL     22644
    ## 4     DUP      1703
    ## 5     DEL      1789
    ## 6     DEL     49256

These dataframe can be provided to the extractSigs\* functions using the
argument `df`. For example:

``` r
extractSigsSv(df=sv_dataframe, vcf.filter='PASS', output='contexts', sv.caller='gridss')
```

### Available signature profiles

Below is a summary of the signature profiles pre-loaded within
mutSigExtractor

  - `SBS_SIGNATURE_PROFILES_V2`: Original 30 SNV signatures
  - `SBS_SIGNATURE_PROFILES_V3`: PCAWG SNV signatures
  - `INDEL_SIGNATURE_PROFILES`: PCAWG indel signatures
  - `DBS_SIGNATURE_PROFILES`: PCAWG DBS signatures
  - `SV_SIGNATURE_PROFILES`: SV signatures from [Nik-Zainal et
    al. 2016](https://www.nature.com/articles/nature17676) with
    clustered SV information removed

### Using other reference genomes

A different reference genome than the default
(BSgenome.Hsapiens.UCSC.hg19) can be used. Genomes should be BSgenomes.
The **variable name** (i.e. no quotes) of the BSgenome object is
specified to `ref.genome`.

``` r
## Make sure to install and load the desired ref genome first
install.packages('BiocManager')
BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')

## Non-default genomes need to be explicitly loaded. The default (BSgenome.Hsapiens.UCSC.hg19)
## is automatically loaded.
library(BSgenome.Hsapiens.UCSC.hg38)

## Specify the name of the BSgenome object to the ref.genome argument
extractSigsSnv(
  vcf.file='/path/to/vcf/with/snvs_and_indels/', vcf.filter='PASS', output='contexts',
  ref.genome=BSgenome.Hsapiens.UCSC.hg38
)
```
