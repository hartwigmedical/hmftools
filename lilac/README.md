# LILAC

## Introduction

LILAC determines the most likely combination of HLA class I alleles present in a patient, aka HLA typing. 
HLA typing is performed to 4 digit allele resolution, meaning LILAC uniquely identifies a specific protein, but ignores synonymous variants 
(6 digits) and intronic differences (8 digits). LILAC is described and validated in the publication:
*Genetic immune escape landscape in primary and metastatic cancer, Nature 2023* ([link](https://www.nature.com/articles/s41588-023-01367-1)).

To get started using LILAC, please jump to [Usage](#usage).

Notable existing tools for HLA class I typing include Polysolver, xHLA, Optitype and DRAGEN-HLA. LILAC offers a number of potential 
advantages including:
- ~99.8% accuracy on 30x-100x WGS samples; ~96% accuracy on panel samples
- Handling of variable sequencing depth
- Reduced false positive detection of rare alleles
- Fully integrated analysis of paired tumor-normal sample data to call allele specific somatic mutations, allelic imbalance, and/or 
complete loss of alleles
- Detection of novel germline variants / alleles (including indels) via analysis of unmatched fragments
- Identification of HLA-Y presence, a pseudogene with high similarity to HLA-A which is present in up to 20% of the population but is not
present in the reference genome

## Table of contents

<!-- TOC -->
* [Introduction](#introduction)
* [Usage](#usage)
  * [Versions](#versions)
  * [Sample inputs](#sample-inputs)
  * [Reference data](#reference-data)
  * [Arguments](#arguments)
* [Output](#output)
  * [Top solution summary](#top-solution-summary)
  * [Top ranked solutions summary](#top-ranked-solutions-summary)
  * [QC metrics](#qc-metrics)
  * [Additional output files](#additional-output-files)
* [Reference data generation](#reference-data-generation)
  * [Allele population frequencies](#allele-population-frequencies)
  * [Allele sequences](#allele-sequences)
* [Algorithm](#algorithm)
  * [Elimination phase](#elimination-phase)
  * [Evidence phase](#evidence-phase)
  * [Tumor and RNA status of alleles](#tumor-and-rna-status-of-alleles)
  * [QC metrics and PON](#qc-metrics-and-pon)
* [Known issues / future improvements](#known-issues--future-improvements)
<!-- TOC -->

## Usage

### Versions

The latest LILAC jar version can be downloaded here: [1.7](https://github.com/hartwigmedical/hmftools/releases/tag/lilac-v1.7.1)

Older versions:
- [1.6](https://github.com/hartwigmedical/hmftools/releases/tag/lilac-v1.6)
- [1.5](https://github.com/hartwigmedical/hmftools/releases/tag/lilac-v1.5.2)
- [1.4](https://github.com/hartwigmedical/hmftools/releases/tag/lilac-v1.4.2)
- [1.3](https://github.com/hartwigmedical/hmftools/releases/tag/lilac-v1.3)
- [1.2](https://github.com/hartwigmedical/hmftools/releases/tag/lilac-v1.2)
- [1.1](https://github.com/hartwigmedical/hmftools/releases/tag/lilac-v1.1)
- [1.0](https://github.com/hartwigmedical/hmftools/releases/tag/lilac-v1.0)

### Sample inputs

LILAC performs HLA typing from GRCh37 or GRCh38 (no alt) BAMs. It can run in various sample modes:
- **Germline-only** or **tumor-only**: Basic HLA typing
- **Paired tumor/normal**: Call allele specific somatic mutations, allelic imbalance, and/or complete loss of alleles. This mode
  can take somatic variant calls and copy number variant calls from [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purple)
  as inputs, but outputs from other callers are also accepted.

LILAC only considers reads from the HLA gene regions. You may provide BAMs which have been sliced for the HLA gene regions.

> [!WARNING]
> BAMs aligned to a GRCh38 genome build containing HLA alt contigs must be realigned to a GRCh38 genome build without HLA alt contigs.
> [AltContigRemapper](https://github.com/hartwigmedical/hmftools/tree/master/bam-tools#altcontigremapper) from hmftools/bam-tools can be
> used for this task.

LILAC has been tested on WGS samples (30-40x germline depth, 100x tumor depth, paired 151 bp reads) and panel samples 
(1000x coverage, paired 86 bp reads). Generally, shorter read length and lower depths are problematic for LILAC. In tumor samples with high purity and LOH, the lost allele in
the tumor may also be difficult to detect.

### Reference data

The reference data required to run LILAC can be downloaded from the `oncoanalyser`
[downloads](https://nf-co.re/oncoanalyser/docs/usage/#reference-data-urls) page:
- Reference genome FASTA
- Allele population frequencies: `lilac_allele_frequencies.csv`
- Allele nucleotide sequences: `hla_ref_nucleotide_sequences.csv`
- Allele amino acid sequences: `hla_ref_aminoacid_sequences.csv`

Allele sequences and frequencies files are found in the `hmf_pipeline_resources.*.tar.gz` bundle from `oncoanalyser`.

See section [Reference data generation](#reference-data-generation) for details on how the allele reference data files are derived.

### Arguments

#### Example commands

##### Reference-only mode

HLA-typing from a reference sample (e.g. blood or normal tissue) BAM.

```
java -jar lilac.jar \
    -sample COLO829R 
    -ref_genome /path/to/ref_genome.fasta \
    -ref_genome_version V37 \
    -resource_dir /path/to/lilac_resources_dir/ \
    -output_dir /output_dir/ \
    -reference_bam COLO829R.bam
```

##### Tumor-only mode

HLA-typing from the tumor sample BAM.

```
java -jar lilac.jar \
    -sample COLO829T \
    -ref_genome /path/to/ref_genome.fasta \
    -ref_genome_version V37 \
    -resource_dir /path/to/lilac_resources_dir/ \
    -output_dir /output_dir/ \
    -tumor_bam COLO829T.bam
```

##### Paired tumor/normal mode

HLA-typing, and calling allele specific somatic mutations, allelic imbalance, and/or complete loss of alleles in the tumor sample. Provide 
the reference and tumor BAMs, as well as somatic variants VCF and gene copy number TSV from 
[PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purple) (or other variant callers).

```
java -jar lilac.jar \
    -sample COLO829T 
    -ref_genome /path/to/ref_genome.fasta \
    -ref_genome_version V37 \
    -resource_dir /path/to/lilac_resources_dir/ \ 
    -output_dir /output_dir/ \
    -reference_bam COLO829R.bam \
    -tumor_bam COLO829T.bam \
    -somatic_vcf COLO829T.purple.somatic.vcf.gz \
    -gene_copy_number COLO829T.purple.cnv.gene.tsv
```

#### Mandatory input paths

| Argument              | Description                                                                                                                                                         |
|:----------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `-sample`             | Sample ID                                                                                                                                                           |
| `-ref_genome`         | Reference genome fasta file                                                                                                                                         |
| `-ref_genome_version` | V37 (default) or V38                                                                                                                                                |
| `-reference_bam`      | Sample germline BAM                                                                                                                                                 |
| `-resource_dir`       | Path to LILAC resource directory containing: <br/>- `hla_ref_aminoacid_sequences.csv`<br/>- `hla_ref_nucleotide_sequences.csv`<br/>- `lilac_allele_frequencies.csv` |

Notes:
- LILAC only considers reads from the HLA gene regions. You may provide BAMs which have been sliced for the HLA gene regions.
- If a tumor BAM is provided in place of the reference BAM, then Lilac will determine the allele solution from the tumor instead.

#### Optional input paths

| Argument            | Description                                                          |
|:--------------------|:---------------------------------------------------------------------|
| `-tumor_bam`        | Sample tumor BAM                                                     |
| `-rna_bam`          | Sample RNA BAM if available                                          |
| `-gene_copy_number` | Sample gene copy number file from PURPLE                             |
| `-somatic_vcf`      | Sample somatic variant VCF file, for annotation of HLA gene variants |

#### Optional parameters

| Argument                          | Default            | Description                                                                                                                                     |
|:----------------------------------|:-------------------|:------------------------------------------------------------------------------------------------------------------------------------------------|
| `-min_base_qual`                  | 30                 | Min base quality for BAM reads, see documentation for details                                                                                   |
| `-min_evidence_support`           | 1                  | Min absolute fragment support for an allele to proceed to evidence phase                                                                        |
| `-min_evidence_factor`            | 0.006              | Min relative fragment support for an allele to proceed to evidence phase, as fraction of total fragments                                        |
| `-min_high_qual_evidence_factor`  | 0.003              | Min relative high base-quality fragment support for an allele to proceed to evidence phase, as fraction of total fragments                      |
| `-min_fragments_per_allele`       | 7                  | See documentation for details                                                                                                                   |
| `-min_fragments_to_remove_single` | 40                 | See documentation for details                                                                                                                   |
| `-top_score_threshold`            | 5                  | Max difference in candidate solution score vs top score as a percentage of total fragments, used for selecting final set of candidate solutions |
| `-hla_y_threshold`                | 0.003              | Min fraction of total reads supporting HLA-Y to consider it present. If present, HLA-Y supporting reads are discarded from HLA typing           |
| `-freq_score_penalty`             | 0.009              | Weight for penalising the score of solutions containing rare alleles                                                                            |
| `-write_types`                    | SUMMARY            | List of output file types to write, separated by `,`. Can be: ALL, SUMMARY, FRAGMENTS, READS, REF_COUNTS                                        |
| `-debug_phasing`                  | Off                | Logs phasing evidence construction                                                                                                              |
| `-expected_alleles`               | Not applied        | List of alleles separated by `;`. These alleles will have their coverage and ranking reported even if not in the winning solution               |
| `-restricted_alleles`             | Not applied        | List of alleles separated by `;`. Restrict evaluation to only these alleles.                                                                    |
| `-threads`                        | 1                  | Number of threads to use for complex evaluation                                                                                                 |
| `-log_level`                      | INFO               | Logs level. Can be: TRACE, DEBUG, INFO, WARN, ERROR                                                                                             |
| `-log_debug`                      | Off (logs at INFO) | Alias for `-log_level DEBUG`                                                                                                                    |

## Output

### Top solution summary

The `<sample_id>.lilac.tsv` file contains a summary of top solution. 
Tumor and/or RNA results are only shown if [tumor and/or RNA inputs](#tumor-and-rna-status-of-alleles) are provided.

| Field(s)                                | Description                                                                                                               |
|:----------------------------------------|:--------------------------------------------------------------------------------------------------------------------------|
| `Allele`                                | Allele ID                                                                                                                 |
| `RefTotal`, `TumorTotal`, `RnaTotal`    | Total assigned fragments from reference/tumor/RNA BAM                                                                     |
| `RefUnique`, `TumorUnique`, `RnaUnique` | Number of fragments assigned uniquely to allele                                                                           |
| `RefShared`, `TumorShared`, `RnaShared` | Number of fragments assigned to allele and others in this solution (i.e. [shared](#shared-fragment-support-across-genes)) |
| `RefWild`, `TumorWild`, `RnaWild`       | Number of fragments matched to wildcard sequences                                                                         |
| `TumorCopyNumber`                       | Copy number from tumor/ref fragment ratio and PURPLE copy number                                                          |
| `SomaticMissense`                       | Matched missense variants                                                                                                 |
| `SomaticNonsenseOrFrameshift`           | Matched nonsense or frameshift variants                                                                                   |
| `SomaticSplice`                         | Matched splice variants                                                                                                   |
| `SomaticSynonymous`                     | Matched synonymous variants                                                                                               |
| `SomaticInframeIndel`                   | Matched inframe indels                                                                                                    |

Example snippet:

```
Allele   RefTotal RefUnique RefShared RefWild TumorTotal ...
A*01:01  6941     4374      2567      0       0          ...
A*24:02  5916     3380      2536      0       0          ...
B*07:02  4103     2222      1881      0       0          ...
B*27:05  3706     1876      1830      0       0          ...
C*02:02  3582     2772      810       0       0          ...
C*07:02  3011     2288      723       0       0          ...
```

### Top ranked solutions summary

The `<sample_id>.lilac.candidates.coverage.tsv` file contains coverage (i.e. fragment support) info for all candidate solutions within X% 
of the top solution's score.

| Field(s)                                    | Description                                                                                                                                           |
|:--------------------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------|
| `Score`                                     | Solution [score](#scoring) (`TotalCoverage` minus penalties)                                                                                          |
| `ComplexityPenalty`, `Complexity`           | [Solution complexity](#scoring) and corresponding penalty                                                                                             |
| `HomozygousCount`                           | Number of homozygous genes                                                                                                                            |
| `CohortFrequencyPenalty`, `CohortFrequency` | [Total log10 population frequency](#scoring) across alleles, and corresponding penalty                                                                |
| `RecoveryCount`                             | Number of alleles in solution that were eliminated but subsequently [recovered](#conditionally-eliminate-and-recover-alleles)                         |
| `WildcardCount`                             | Number of wildcards characters in all allele sequences                                                                                                |
| `TotalCoverage`                             | Sum of total fragments assigned across all alleles in solution                                                                                        |
| `UniqueCoverage`                            | Sum of uniquely assigned fragments across all alleles in solution                                                                                     |
| `SharedCoverage`                            | Sum of [shared](#shared-fragment-support-across-genes) fragments across all alleles  in solution                                                      |
| `WildCoverage`                              | Sum of fragments matched to wildcard sequences across all alleles in solution                                                                         |
| `A1`, `A2`, `B1`, `B2`, etc...              | Per gene (e.g. HLA-`A`) and allele slot (e.g. `1`/`2`), allele name and `[total, unique, shared, wildcard]` fragment counts per allele                |

Example snippet:

```
Score    ComplexityPenalty Complexity HomozygousCount CohortFrequencyPenalty CohortFrequency
15900.63 -1060.57          31         0               -144.80                -9.41          
15884.65 -1060.01          31         0               -152.34                -9.90          
...      ...               ...        ...             ...                    ...

RecoveryCount WildcardCount TotalCoverage UniqueCoverage SharedCoverage WildCoverage
0             0             17106         8171           8935           0           
1             0             17097         8196           8901           0           
...           ...           ...           ...            ...            ...

A1                        A2                        B1                       ...
A*02:01[3642,2326,1316,0] A*33:03[4332,3034,1298,0] B*51:01[2067,748,1319,0] ...
A*02:01[3642,2326,1316,0] A*33:03[4332,3034,1298,0] B*51:12[2041,739,1302,0] ...
...                       ...                       ...                      ...
```

### QC metrics

The `<sample_id>.lilac.qc.tsv` file contains comprehensive QC metrics.

| Field(s)                                                         | Description                                                                                                                                                                             |
|:-----------------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `Status`                                                         | Either PASS or 1 or more warnings (see below for warning descriptions)                                                                                                                  |
| `ScoreMargin`                                                    | Difference in score to second-top solution                                                                                                                                              |
| `NextSolutionAlleles`                                            | Allele difference in second-top solution                                                                                                                                                |
| `MedianBaseQuality`                                              | Median base quality across all coding bases from all fragments                                                                                                                          |
| `HlaYAllele`                                                     | HLA-Y allele detected if present, else 'NONE'                                                                                                                                           |
| `DiscardedIndels`                                                | Discarded fragments due to unknown indels                                                                                                                                               |
| `DiscardedIndelMaxFrags`                                         | Maximum fragment support for an indel detected but not present in any known allele.                                                                                                     |
| `DiscardedAlignmentFragments`                                    | Fragments discarded because 1 read aligns more than 1000 bases from an HLA gene                                                                                                         |
| `A_LowCoverageBases`, `B_LowCoverageBases`, `C_LowCoverageBases` | Number of bases with less than 15-depth coverage across all coding bases                                                                                                                |
| `ATypes`, `BTypes`, `CTypes`                                     | Number of distinct alleles fitted (0,1 or 2) per HLA gene                                                                                                                               |
| `TotalFragments`                                                 | Total of fragments overlapping a coding base of HLA-A, HLA-B or HLA-C with MAPQ >=1                                                                                                     |
| `FittedFragments`                                                | Total fragments assigned to fitted alleles (ie. exact or wildcard match at every heterozygous location)                                                                                 |
| `UnmatchedFragments`                                             | Fragment that do not match any allele exactly in solution at all heterozygous location (eg. due to sequencing error or mapping error)                                                   |
| `UninformativeFragments`                                         | Fragments that do not overlap any heterozygous location considered in evidence phase                                                                                                    |
| `HlaYFragments`                                                  | Fragments excluded from fit due to allocation to assignment to HLA-Y.                                                                                                                   |
| `PercentUnique`                                                  | Percentage of fitted fragments that are uniquely assigned to 1 allele                                                                                                                   |
| `PercentShared`                                                  | Percentage of fitted fragments allocated across multiple alleles                                                                                                                        |
| `PercentWildcard`                                                | Percentage of fitted fragments uniquely assigned to wildcard regions                                                                                                                    |
| `UnusedAminoAcids`                                               | # of amino acids detected with at least 3 fragments support but not present in final allele solution set (excluding PON filtered and regions overlapping wildcards in selected alleles) |
| `UnusedAminoAcidMaxFrags`                                        | Maximum fragments support for an unmatched amino acid not present in the final allele solution set                                                                                      |
| `UnusedHaplotypes`                                               | # of haplotypes phased with at least 3 fragments support but not present in final allele solution set (excluding PON filtered and regions overlapping wildcards in selected alleles)    |
| `UnusedHaplotypeMaxFrags`                                        | Maximum support for an unmatched haplotype not present in final allele solution set                                                                                                     |
| `SomaticVariantsMatched`                                         | Somatic variants supported by solution allele                                                                                                                                           |
| `SomaticVariantsUnmatched`                                       | Somatic variants not supported by solution allele                                                                                                                                       |

Example snippet:

```
Status ScoreMargin NextSolutionAlleles MedianBaseQuality HlaYAllele ...
PASS   26.513      C*02:175            37                NONE       ...
```

### Additional output files

| File                                           | Description                                                                                                       |
|:-----------------------------------------------|:------------------------------------------------------------------------------------------------------------------|
| `<sample_id>.lilac.log`                        | Log file with information about the fit and details of all [unmatched haplotypes and indels](#qc-metrics-and-pon) |
| `<sample_id>.lilac.HLA-A.aminoacids.txt`       | Fragment support for each amino acid by HLA gene                                                                  |
| `<sample_id>.lilac.HLA-A.nucleotides.txt`      | Fragment support for each nucleotide by HLA gene                                                                  |
| `<sample_id>.lilac.fragments.csv`              | Read details for all BAM fragments                                                                                |
| `<sample_id>.lilac.candidates.fragments.csv`   | Allocation of each fragment to one or more solutions and which alleles they support                               |
| `<sample_id>.lilac.candidates.aminoacids.txt`  | Fragment support for amino acids in the candidate alleles                                                         |
| `<sample_id>.lilac.candidates.nucleotides.txt` | Fragment support for nucleotides in the candidate alleles                                                         |
| `<sample_id>.lilac.somatic.vcf.gz`             | Annotation of supplied somatic variants if somatic VCF provided with assigned allele                              |

## Reference data generation

### Allele population frequencies

#### Web scraping of AFND

Allele population frequencies were scraped from the Allele Frequency Net Database (AFND; http://www.allelefrequencies.net/) on 2025-01-28
using the [hla_allele_freq_downloader.py](./src/main/resources/frequencies/hla_allele_freq_downloader.py) Python script:

```
python hla_allele_freq_downloader.py \
  --output_dir /output_dir/ \
  --locus A,B,C \
  --hla_level 4 \
  --hla_locus_type Classical \
  --population_standard g
```

#### Normalising frequencies

The AFND entry represents the frequency of an allele in a population (typically country or ethnicity). However, each frequency may have been 
calculated as part of a larger study. Since sample sizes and method of frequency calculation vary between studies, we calculated normalised 
frequencies using the [hla_allele_freq_normaliser.py](src/main/resources/frequencies/hla_allele_freq_normaliser.py) Python script:

```
python hla_allele_freq_normaliser.py \
--input_file afnd.hla.gold_standard.tsv.gz \
--output_dir /output_dir/ \
--write_debug_files
```

This script performs the below procedure:
- Exclude entries where frequency is 0 
- Calculate approx_observations (= sample_size * raw_allele_frequency), and exclude entries where this is ≤2
- Calculate the mean frequency **per allele**, weighted by study cohort size, where study cohort size is capped at 1000 samples
- Scale frequencies to add to 1 **per HLA gene**

### Allele sequences
Allele sequences were derived from the IMGT/HLA database (downloaded on 2024-12-19), accessible with these links:
- https://www.ebi.ac.uk/ipd/imgt/hla/download/ or
- https://github.com/ANHIG/IMGTHLA

LILAC constructs its own sequence files using the nucleotide multiple sequence alignment files from the database:
`A_nuc.txt`, `B_nuc.txt`, `C_nuc.txt`, `Y_nuc.txt`. The procedure is described below.

#### Wildcard handling

> [!NOTE]
> This is currently performed using a Python script; to be implemented in Java in `GenerateReferenceSequences`

Wildcards (`*`) are present in some allele nucleotide sequences particularly in exon 1 and exons 4-8 (e.g. due to incomplete sequencing). 
The nucleotide sequence of wildcard sequences are inferred by first generating a consensus sequence per 2-digit allele group.

First we determine the representative set of 4-digit alleles. If a 4-digit allele is only represented by its 6-digit (synonymous variants)
or 8-digit (intronic variants) alleles we greedily select the first (e.g. for `A*34:01` we select `A*34:01:01:01` from 
`A*34:01:01:01`, `A*34:01:01:02`, `A*34:01:02`, ...)

For each 2-digit allele group, (e.g. group `A*34` is composed of alleles `A*34:01`, `A*34:02`, `A*34:03`, ...), we select the 4-digit 
alleles fulfilling the below criteria:
- ≥0.001 population frequency
- Has no wildcards
- Has no expression suffix (last letter in allele name, e.g. `C*04:09N`)

As a fallback, if no 4-digit alleles fulfill the above criteria, we select the most common. This usually occurs when all alleles in the 
group have wildcards. Note, in case of this fallback condition, the most common allele is by definition the consensus sequence.

A consensus sequence is then created from the selected 4-digit allele sequences, with wildcards being are assigned if alleles in the 
2-digit group have conflicting nucleotides. For example, for the first few nucleotides of `A*34`:

```
A*34:01 sequence: ATG GCC ATC
A*34:02 sequence: ATG GCC GTC
A*34 consensus:   ATG GCC *TC
```

For every allele nucleotide sequence, wildcard nucleotides are then replaced with the nucleotide from the corresponding 2-digit group 
consensus sequence, and these new sequences written to new `*_nuc.txt` files. For example, for the first few nucleotides of `A*34:01:02`:

```
A*34:01:02 original: *** *** ***
A*34:01:02 replaced: ATG GCC *TC
```

#### Generate LILAC sequences

These files are converted to the LILAC allele sequences files with the following command:
```
java -cp lilac.jar com.hartwig.hmftools.lilac.utils.GenerateReferenceSequences \
   -resource_dir /dir_containing/*.txt \ 
   -output_dir /output_dir/
```

This routine first selects the representative set of 4-digit alleles. If a 4-digit allele is only represented by its 6-digit 
(synonymous variants) or 8-digit (intronic variants) alleles we greedily select the first (e.g. for `A*34:01` we select `A*34:01:01:01` 
from `A*34:01:01:01`,  `A*34:01:01:02`, `A*34:01:02`, ...).

Then, the nucleotide multiple sequence alignments are converted into nucleotide and amino acid sequences, 
with deletions represented by a dot (`.`). For example, for the beginning sequence of `A*30:220Q`:

```
Reference sequence:    ATG GCC GTC ATG GCG CCC CGA ACC CTC CTC CTG CTA CTC
A*30:220Q alignment:   --- --- --- --- --- --- --- --- --. ... ... ..- ---

# Spaces shown for here clarity but absent in the LILAC nucleotides file
A*30:220Q nucleotides: ATG GCC GTC ATG GCG CCC CGA ACC CT. ... ... ..A CTC 

A*30:220Q amino acids: MAVMAPRTLL...SG
```

## Algorithm

The starting data for the LILAC algorithm is:
- HLA-A, HLA-B, HLA-C and HLA-Y [4-digit allele sequences](#allele-sequences), excluding the following due to their high similarity closely 
related pseudogenes:
  - HLA-H: `A*31:135`, `A*33:191`, `A*02:783`, `B*07:282`
  - HLA-Y: `A*30:205`, `A*30:207`, `A*30:225`, `A*30:228`, `A*01:81`, `A*01:237`
- All fragments in a BAM aligned to HLA-A, HLA-B and HLA-C gene regions, which:
  - Are not duplicates
  - Have at least 1 read with an alignment overlapping a coding base of HLA-A, HLA-B or HLA-C
  - Have all alignments within 1000 bases of an HLA coding region
  - Have a mapping quality of at least 1

The main steps of the algorithm are:
- Elimination phase - eliminate candidate alleles: 
  - Construct **nucleotide** matrix (position x fragment support per base). Eliminate alleles with sequences not matching possible **nucleotide** sequences
  - Construct **amino acid** matrix (position x fragment support per residue). Eliminate alleles with sequences not matching possible **amino acid** sequences
  - Eliminate based on **phased haplotypes**
  - Exclude HLA-Y pseudogene fragments
  - Conditionally eliminate and recover alleles
- Evidence phase - score and rank allele combinations

After the germline alleles are determined, LILAC determines the tumor copy number and any somatic mutations in each allele. Note that if
more than 300 bases of the HLA-A, HLA-B and HLA-C coding regions have less than 10 coverage, then LILAC will fail with errors and will not
try to fit the sample.

### Elimination phase
The elimination phase removes alleles that are unlikely part of the final solution, namely, if they do not have sufficient fragment support 
at each nucleotide / amino acid position. This reduces runtime by reducing the number of allele combinations LILAC needs to consider in 
the evidence phase.

#### Nucleotide matrix construction

For each HLA gene, create a matrix containing the nucleotide counts (i.e. high quality fragment support) for each position, resulting in the
below output. Bases with 0 fragment support not shown, and soft-clipped bases are represented as multi-character strings.

```
#index count1 base1 count2 base2 etc...
0	766	A	991
1	769	T	997
2	768	G	999
3	774	G	772	C	237
4	776	C	773	G	167	T	74
5	778	C	776	G	237
6	786	G	1025
7	792	T	1034
8	793	C	1033
9	803	A	1048
...
286	1167	A	3682
287	1181	G	3120	A	522	AATAAGCAAAACAAACACACAA	55
288	1179	A	3412	G	289
...
```

Positions with more than 1 amino acid candidate are deemed **heterozygous**, and positions with only 1 are considered to be **homozygous**.

Construction of the nucleotide matrix has several steps / conditions which are described in the below subsections.

##### Shared fragment support across genes

Firstly, fragment bases aligned at exons with boundaries identical across two or more HLA genes will count towards fragment support for the 
corresponding nucleotide positions and genes. For example, the below table shows the exon boundaries (nucleotide end positions) for each 
HLA class I gene. Fragment bases at exons 1-4 will contribute to fragment counts in the HLA-A, HLA-B and HLA-C matrices, whereas fragment 
bases at exons 5-6 will support only contribute to the HLA-A/B matrices or the HLA-C matrix.

| Gene \ Exon | 1  | 2   | 3   | 4   | 5    | 6    | 7    |
|-------------|----|-----|-----|-----|------|------|------|
| HLA-A       | 72 | 342 | 618 | 894 | 1011 | 1044 | 1092 |
| HLA-B       | 72 | 342 | 618 | 894 | 1011 | 1044 | 1086 |
| HLA-C       | 72 | 342 | 618 | 894 | 1014 | 1047 | 1095 |

##### Indel handling

Fragments with in-frame indels only contribute to counts in the nucleotide matrix if the indel matches an existing HLA allele allowing for 
realignment. Fragments with out-of-frame indels are always excluded.

##### Filtering for nucleotide candidates

For each position, candidate nucleotide bases are filtered for those with fragment support greater than or equal to:
- `min_hi_qual_fragment_support` = `max(1, min_high_qual_evidence_factor * hi_qual_fragment_support)` where:
  - `min_high_qual_evidence_factor`: 0.003 by default
  - `hi_qual_fragment_support`: number of fragment bases with at least `min(30, median_base_quality)` supporting a candidate base

- `min_overall_fragment_support` = `max(2, min_evidence_factor * overall_fragment_support)` where:
  - `min_evidence_factor`: 0.006 by default
  - `overall_fragment_support`: total number of fragment bases supporting a candidate base

#### Elimination based on nucleotide matrix

Any alleles with a nucleotide or inframe indel that do not match the possible sequences of nucleotide matrix are eliminated.

For example, `A*01:237` would be eliminated because the `G` at index 9 does not match the `A` in the above example 
[nucleotide matrix](#nucleotide-matrix):

```
Sequence: ATGGCCGTCG...
Index:    0123456789...
```

An exception is that if a position in the nucleotide matrix has <10 total fragment support, that position will not be used to eliminate 
alleles. This prevents inadvertently eliminating alleles due to low depth as a result of sequencing issues, e.g. due certain SNVs reducing 
probe binding affinity.

#### Amino acid matrix

Similar to the nucleotide matrix, LILAC also constructs a matrix of amino acid candidates. Amino acids with 0 fragment support not shown, 
and amino acids derived from soft-clipped bases are represented as multi-character strings.

```
0	749	M	971
1	755	A	752	R	164	L	72
2	776	V	1009
3	792	M	1019	T	14
...
94	1103	S	1542	T	1090	A	936
95	1155	Q	3556	QISKTNTQ	55
96	1178	T	3385	A	285
...
```

Positions with more than 1 amino acid candidate are deemed **heterozygous**, and positions with only 1 are considered to be **homozygous**.

The same steps used in nucleotide matrix construction are used in amino acid matrix construction, namely: 
- [Calculating shared fragment support across genes](#shared-fragment-support-across-genes)
- [Indel handling](#indel-handling)
- [Filtering of amino acid candidates](#filtering-for-nucleotide-candidates)

An addition step, exon boundary 'enrichment', is applied to exons 1-4 where [exon boundaries are shared](#shared-fragment-support-across-genes)
(i.e. before nucleotide index 894, amino acid index 298). For codons crossing these exon boundaries, any associated fragment is enriched 
with **homozygous** nucleotide candidates on the other side of the exon boundary so that an amino acid can constructed.


#### Elimination based on amino acid matrix

Alleles are eliminated based on the amino acid matrix in the same way as for the 
[nucleotide matrix](#elimination-based-on-nucleotide-matrix), i.e. any alleles with an amino acid or inframe indel that do not match the 
possible sequences of amino acid matrix are eliminated.

#### Elimination based on phased haplotypes

In this step, we identify 'phased haplotypes', consecutive **heterozygous** amino acids ('haplotype') supported by contiguously
overlapping fragments ('phased'). We eliminate any alleles whose (sub)sequences do not match with any identified phased haplotype.

First, we find phased evidence of each consecutive **pair** of heterozygous amino acids and record the haplotypes of all fragments 
containing both amino acids:
- This is performed separately for each of HLA-A, HLA-B and HLA-C genes to account for differences in exon 
boundaries. Fragments overlapping amino acid index 298 (i.e. [nucleotide index 894](#shared-fragment-support-across-genes)) onwards will 
only be assignable to a subset of the alleles, since the exon boundaries differ after this amino acid.
- The points are only phased together if the total coverage is at least 7 fragments per allele (`min_fragments_per_allele`) included in the 
subset at that location . This effectively means between 14 fragments (for 2 alleles) and fragments 42 (for 6 alleles) with shared 
coverage depending on the amino acid location.
- A phased haplotype with only 1 supporting fragment will be removed if the total fragments supporting the pair is 
40 or more (`min_fragments_to_remove_single`) as it is assumed to be a sequencing error.

We then iteratively choose the phased haplotype with the most support and perform the following routine:

1) Find other phased evidence that overlaps it.

2) Find the minimum number of amino acid locations required to uniquely identify each phased evidence. For example, if the left evidence has 
haplotypes [SP, ST] you only need the last amino acid but if the right evidence has haplotypes [PD, TD, TS] you would need both amino acid 
locations.

3) Find evidence of fragments that contain all the required amino acids. As with the paired evidence, there must be at least 7 fragments per 
allele supporting the pair (`min_fragments_per_allele`) and a haplotype with only 1 supporting fragment will be removed if the total 
fragments supporting the pair is 40 or more (`min_fragments_to_remove_single`).

4) Check that the new overlapping evidence is consistent with the existing evidence.

5) Merge the new evidence with the existing paired evidence.

6) Replace the two pieces of used evidence with the new merged evidence.

Once complete, we eliminate any alleles that do not match the phased evidence. 

#### Excluding HLA-Y pseudogene fragments
HLA-Y is a pseudogene that is highly similar to HLA-A, but is not present in the human ref genome, and is found in approximately 17% of the
Hartwig cohort. The presence of HLA-Y can cause confusion in HLA typing, particularly of HLA-A alleles.

HLA-Y is considered present in a sample if **>0.3%** of fragments align uniquely to any of the HLA-Y allele sequences. If HLA-Y is present,
**any** fragment which matches exactly to an HLA-Y allele (uniquely or shared with other alleles) are excluded from further analysis.

#### Conditionally eliminate and recover alleles

The below steps are performed to further eliminate alleles, while preventing inadvertently eliminating a potential true allele.

1) Identify high confidence 2-digit allele groups (**>2%** unique fragment support). Discard any 4-digit alleles not belonging to the high 
confidence 2-digit allele groups.

2) Identify low confidence 2-digit allele groups (**>0.1%** unique fragment support, or >3.5% of total fragment support). Recover maximally 
two common 4-digit alleles per HLA gene from each group. A common allele is defined as having 
(>0.0001 [population frequency](#normalising-frequencies)). 

3) For each remaining 4-digit allele that is 1) not common or 2) has wildcards (`*`) in its sequence, determine the 2-digit allele group and 
recover common 4-digit alleles belonging to this group.

4) Recover `C*04:09N`, the most common HLA allele with a frameshift variant leading to a stop loss

5) Identify high confidence 4-digit alleles (**>1%** unique fragment support). These alleles **must** be a part of the allele combinations
considered in the evidence phase.


### Evidence phase

From the remaining 4-digit alleles, all 'complexes' (i.e. candidate allele combinations) are evaluated, with each 
complex requiring:
- 1 or 2 alleles assigned per gene (as each gene can be homozygous or heterozygous)
- At least 1 allele must match each uniquely supported 2-digit allele group

#### Determining fragment support per complex

For each complex, LILAC counts the number of fragments that can be aligned exactly to at least one allele in the complex at all heterozygous 
locations (i.e. where the amino acids differ between alleles).
- When checking if a fragment can be aligned to an allele, if a fragment has an amino acid which does not match **any** of the amino acid
candidates at a heterozygous position, but has at least 1 nucleotide with base quality >`min(30, median_base_quality)`, then the amino 
acid is deemed to match all amino acid candidates.
- Any fragments that can be aligned to multiple alleles are apportioned equally between the alleles and counted as shared fragments.
- For exon boundaries only exact nucleotide matches are permitted. 

Wildcards in allele sequences are handled as below:
- Fragments which don't match the exact (sub)sequence of any candidate allele are dropped altogether from consideration
- For each remaining fragments, a fragment is considered to support an allele if the fragment sequence matches all non-wildcard sequences 
of that allele (i.e. any amino acid is deemed a match to a wildcard)

#### Scoring

The score for a complex is calculated based on `total_coverage` (sum of the fragments aligned to each allele), and several penalties: 

```
score = total_coverage - frequency_penalty - solution_complexity_penalty - recovery_penalty
```

`frequency_penalty` penalises solutions with rare alleles:
```
frequencies = [0.112, 0.019, NaN]
frequencies_filled = [0.112, 0.019, 0.0001] # min_population_frequency applied to alleles without frequency data
sum_log_frequency = sum(-log10(frequencies_filled))

frequency_penalty = total_coverage * sum_log_frequency * 0.009 # frequency_penalty_weight
```

`solution_complexity_penalty` penalises solutions with more alleles (i.e. favor simpler solutions):

```
# Allele sequences in complex:
| exon    | 1    | 2    | 3    | 4    | 5    | 6    | 7    |
| allele1 | MAVM | SHSM | SHTV | APKT | LSSQ | RKGG | SDSA |
| allele2 | MRVT | SHSM | SHII | PPKT | LSSQ | GKGG | SDSA |
| allele3 | MRVM | SHSM | SHTL | HPKT | LSSQ | GKGG | SNSA |

unique_sequences_per_exon = [ 3, 1, 3, 3, 1, 2, 2 ]
solution_complexity = sum(unique_sequences_per_exon)

solution_complexity_penalty = total_coverage * solution_complexity * 0.002 # solution_complexity_penalty_weight
```

`recovery_penalty` penalises solutions having recovered alleles (note: `recovery_penalty` is currently turned off):

```
recovery_penalty = total_coverage * recovered_alleles_count * 0 # recovery_penalty_weight
```

Two performance optimisations are made at this stage of scoring:
- If there are predicted to be >1,000,000 complexes, the evidence phase is first performed individually for complexes of 2
alleles per gene to find the top candidates for each gene. LILAC retains only the top 5 pairs including each individual allele candidate and 
then chooses the first 10 unique alleles appearing in the ranked list of pairs, with any common alleles also retained. The evidence phase 
is then subsequently run using this reduced set of candidate alleles.
- Complexes are scored with the number of fragments is downsampled to a maximum of 10,000

Then, complexes within `top_score * 0.005` (`top_score_threshold`) are rescored with all fragments (i.e. without downsampling) to calculate 
the actual complex scores and rankings. 

The complex with the top score is chosen as the final solution. If 2 or more complexes have the exact same score, the solution with the 
lowest alphabetical 4-digit alleles is chosen.

### Tumor and RNA status of alleles

#### Somatic variant assignment to alleles

LILAC optionally accepts a somatic SNV/indel VCF (e.g. from 
[SAGE](https://github.com/hartwigmedical/hmftools/tree/master/sage) or 
[PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purple)) 
to assign somatic variants to the specific allele which is damaged.

LILAC gathers the set of PASS variants from the VCF that overlap either a coding or canonical splice region in any of HLA-A,
HLA-B and HLA-C and finds all fragments that contain that variant. LILAC assigns apportions the fragment to the allele which matches the
fragment at all heterozygous locations after excluding any somatic variants from the fragment. The allele with the highest matching
fragment count is determined to contain the somatic variant. If the variant is assigned to 2 alleles with identical weight, it is
assigned with 0.5 weight to each. In the case of homozygous alleles, the variant is assigned arbitrarily to the 1st allele.


#### Tumor allele specific copy number

LILAC optionally accepts a tumor BAM file and a [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purple) 
gene copy number file (containing minimum copy number and minimum minor allele copy number per gene) to calculate tumor allele specific 
copy number.

If these inputs are provided, then fragments are counted for each allele in the determined type. For each HLA gene (HLA-A,
HLA-B, HLA-C), the minor allele copy number is assigned to the allele with the lowest ratio of supporting fragments for the allele in the 
tumor compared to the normal. The other allele for each is assigned the implied major allele copy number from the gene copy number file. 
If a gene is homozygous present in the normal sample, then the minor and major allele copy numbers are arbitrarily assigned in the tumor.

#### RNA 'expression' of alleles

LILAC also optionally accepts an RNA BAM. Similar to determining [tumor allele specific copy number](#tumor-allele-specific-copy-number) 
fragments in the BAM are counted for each allele in the determined type. This can be interpreted as a proxy for allele specific expression.

### QC metrics and PON

LILAC produces a comprehensive set of QC metrics and provides warning statuses if the typing confidence may be diminished. In general, the 
warnings may indicate the presence of a novel germline variant/allele, presence of an unusual pseudogene type, or an incorrectly typed 
sample.

The table below shows the possible QC metrics:

| Warning                        | Description                                                                                                          |
|:-------------------------------|:---------------------------------------------------------------------------------------------------------------------|
| WARN_UNMATCHED_TYPE            | all A, B or C types eliminated                                                                                       |
| WARN_UNMATCHED_INDEL           | A novel (non PON) indel is present that was not fit to any allele supported by at least 0.5% of fitted fragments.    |
| WARN_UNMATCHED_HAPLOTYPE       | A novel (non PON) haplotype is present that was not fit to any allele supported by at least 1% of fitted fragments.  |
| WARN_UNMATCHED_AMINO_ACID      | A novel (non PON) amino acid is present that was not fit to any allele supported by at least 1% of fitted fragments. |
| WARN_UNMATCHED_SOMATIC_VARIANT | A somatic variant from the input vcf could not be assigned to any allele                                             |
| WARN_LOW_BASE_QUALITY          | Median base quality < 25                                                                                             |
| WARN_LOW_COVERAGE              | More than 50 distinct bases of HLA-A, HLA-B and HLA-C combined have less than 10 coverage.                           |
| FAIL_LOW_COVERAGE              | More than 200 distinct bases of HLA-A, HLA-B and HLA-C combined have less than 10 coverage.                          |

In the case of unmatched haplotypes and indels, we find that there are a number of recurrent artefacts found across many samples in our 
cohort. We therefore exclude the following indels and haplotypes prior to calculating QC metrics for the fit.

#### Indel PON

| 37_Position | 38_Position | Variant | Curation                 |
|:------------|:------------|:--------|:-------------------------|
| 29912028    | 29944251    | AN>A    | Recurrent HLA-H artefact |
| 29911899    | 29944122    | A>ACC   | Recurrent HLA-Y artefact |
| 29911318    | 29943541    | CN>C    | Recurrent artefact       |
| 31324524    | 31356747    | GNN>G   | Phased alignment issue   |
| 31324528    | 31356751    | G>GGT   | Phased alignment issue   |
| 29910716    | 29942939    | CN>C    | Recurrent artefact       |
| 29910657    | 29942880    | G>GT    | Recurrent artefact       |
| 29912392    | 29944615    | A>AGG   | Recurrent artefact       |
| 29911899    | 29944122    | A>AC    | Recurrent artefact       |

#### Haplotype PON

| Start | Haplotype                                                                                 | Curation                         |
|:------|:------------------------------------------------------------------------------------------|:---------------------------------|
| 1     | AVVAPRTLLLLLSGALALTQTWAG                                                                  | HLA-Y                            |
| 25    | SHSMRYFSTSVSRPGSGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWMEQEEPEYWDRQTEISKTNAQIDLESLRIALRYYNQSED | HLA-Y                            |
| 47    | AVGYVDDTQFVRFDSDAASQRMEPRAPWMEQEEPEYWDRQTQISKTNAQIDLESLRIALR                              | HLA-Y                            |
| 96    | IDLESLRIALRYYNQSEA                                                                        | HLA-Y                            |
| 117   | TIQRMSGCDVGSDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITQRKWEAARQAEQLRAYLEGECMEWLRRYLENGKETLQRT | HLA-Y                            |
| 206   | DAPKTHMTHHAVSDNEATLRCXALSFYPAEITLTWQR                                                     | HLA-Y                            |
| 206   | DPPKTHMTHYPISDHEATLRCWALG                                                                 | HLA-H                            |
| 206   | DPPKTHMTHHPISDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDG                                 | HLA-H                            |
| 206   | DPPKTHMTHHPISDHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGDG                                 | Likely HLA-Y                     |
| 242   | RDGEDQTQDTELVETRPAGDGTFQKWASVVVPSGQEQRYTCHVQHEGLPKPLTLRWE                                 | Likely HLA-Y                     |
| 262   | GIFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWE                                                     | HLA-Y                            |
| 298   | EPSSQPTIPIVGILAGLVLFGAVIAGAVVAAVMWRRKS                                                    | Likely HLA-Y                     |
| 299   | PSSHPTIPIVGILAGLVLFGAVIAGAVVAAVMWRRKS                                                     | HLA-Y                            |
| 299   | IPNLGIVSGPAVLAVLAVLAVLAV                                                                  | Alignment issue with C\*17 indel |
| 337   | DRKGGSYSQAAS                                                                              | HLA-Y or HLA-H                   |
| 351   | IAQGSDVSLTAC                                                                              | HLA-Y                            |

## Known issues / future improvements

- **HG38 reference genome with alt contigs are not supported** - LILAC currently only obtains reads aligned to the HLA Class 1 genes.  Also need to get all the contigs with ref_name =~ /^HLA|chr6.\*alt/"
- **LILAC should obtain all read with mates that map at or near HLA class I genes** - Some reads may have been mismapped to pseudogenes or HLA-H and can be rescued by their mates and remapped to HLA class I genes
- **Support for explicit novel allele prediction** - LILAC reports unmatched_haplotype, but do not explicitly predict the AA sequence of novel alleles.
- **Allele elimination in tumor samples with high purity & LOH** -  Where there is a LOH in the tumor combined with a high purity one of the germline alleles may have very little support in the tumor.  In such cases it is possible that LILAC eliminate the correct allele and end up choosing another similar but incorrect allele (generally either a rare allele which escaped elimination or a common allele that was recovered), instead of calling the LOH.   Depending on the coverage pattern, this  could occur at amino acid elimination or during phasing.   It may make sense to use lower thresholds at the heterozygous amino acid stage and impose higher minimum coverage requirements in the phasing stage to counter this.
- **Recovery** - LILAC currently only recover a common allele if the allele is part of a uniquely supported 2 digit type.   However where alleles are similar across groups, we may not always find unique support for a single 2 digit type.
- **Elimination of alleles where coverage is very weak** - If coverage is very weak for a specific base (say <8 reads) in a gene we should mark that location as heterozygous AND not phase or eliminate alleles based on that location.  This may impact low or variable coverage samples, especially where base quality is low.
- **Indel realignment** - Sometimes we miss evidence supporting INDELs due to realignment issues.  In particular if a fragment does not have a soft clip or an INDEL it is not considered for realignment.  This can be a problem for C\*17 indels which have long homology and sometimes align without soft clip or INDEL.
- **Generic handling of out-of-frame indel** - We currently have special handling for C\*04:09N only.    At least 2 other alleles are relatively common: B\*51:11N, A\*24:09N.   We should handle this more generically.
- **Quality trimming** - We currently quality trim all heterozygous amino acids which have at least 1 nt with base qual <min(30,medianBaseQuality).   A more optimal logic would be to fuzzy match the amino acids.  A better algorithm would be the following:
 <pre>
If a fragment has an amino acid which does not match ANY of the amino acid candidates at a heterozygous location, but has at 
least 1 nucleotide with base qual>=min(30,medianBaseQuality), then the amino acid is deemed to match if it matches any 
combination of the exact high confidence nucleotides with base qual >=min(30,medianBaseQuality) together with any of the 
candidate nucleotides at the bases locations with qual < min(30,medianBaseQuality).   For exon boundaries only exact 
nucleotide matches are permitted.
 </pre>


 
