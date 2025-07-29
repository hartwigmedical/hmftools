# LILAC

## Introduction

LILAC determines the most likely combination of HLA class I alleles present in a patient, aka HLA typing. 
HLA typing is performed to 4 digit allele resolution, meaning LILAC uniquely identifies a specific protein, but ignores synonymous variants 
(6 digits) and intronic differences (8 digits). LILAC is described and validated in the publication:
*Genetic immune escape landscape in primary and metastatic cancer, Nature 2023* ([link](https://www.nature.com/articles/s41588-023-01367-1)).

To get started using LILAC, please [download the reference data](#reference-data), and see section [Usage](#usage) to download the LILAC jar
and see the run commands.

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

## Sample inputs
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

LILAC has been optimised on:
- WGS samples: 30-40x germline depth, 100x tumor depth, paired 151bp reads
- Panel samples: 1000x coverage, paired 86bp reads

Generally, shorter read length and lower depths are problematic for LILAC. In tumor samples with high purity and LOH, the lost allele in 
the tumor may also be difficult to detect.

## Reference data

The reference data required to run LILAC can be downloaded from the `oncoanalyser` 
[downloads](https://nf-co.re/oncoanalyser/docs/usage/#reference-data-urls) page:
- Reference genome FASTA
- Allele population frequencies: `lilac_allele_frequencies.csv`
- Allele nucleotide sequences: `hla_ref_nucleotide_sequences.csv`
- Allele amino acid sequences: `hla_ref_aminoacid_sequences.csv`

Allele sequences and frequencies files are found in the `hmf_pipeline_resources.*.tar.gz` bundle from `oncoanalyser`.

See section [Reference data generation](#reference-data-generation) for details on how the allele reference data files are derived.

## Usage

The LILAC jar version can be downloaded [here](https://github.com/hartwigmedical/hmftools/releases/tag/lilac-v1.7.1). Older LILAC jar 
versions can be downloaded from [hmftools > releases](https://github.com/hartwigmedical/hmftools/releases?q=lilac).

### Examples

#### Reference-only mode

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

#### Tumor-only mode

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

#### Paired tumor/normal mode

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

### Mandatory input paths

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

### Optional input paths

| Argument            | Description                                                          |
|:--------------------|:---------------------------------------------------------------------|
| `-tumor_bam`        | Sample tumor BAM                                                     |
| `-rna_bam`          | Sample RNA BAM if available                                          |
| `-gene_copy_number` | Sample gene copy number file from PURPLE                             |
| `-somatic_vcf`      | Sample somatic variant VCF file, for annotation of HLA gene variants |

### Optional parameters

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


## Reference data generation

### Allele population frequencies

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

Wildcards (`*`) are present in some allele nucleotide sequences (e.g. due to incomplete sequencing). The nucleotide sequence of
wildcard sequences are inferred by first generating a consensus sequence per 2-digit allele group.

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
2-digit group have conflicting nucleotides, for example in `A*34`:

```
A*34:01:        ATG GCC ATC ...
A*34:02:        ATG GCC GTC ...
A*34 consensus: ATG GCC *TC ...
```

For every allele nucleotide sequence, wildcard nucleotides are then replaced with the nucleotide from the corresponding 2-digit group 
consensus sequence, and these new sequences written to new `*_nuc.txt` files. For example, the sequence of `A*34:01:02` becomes:

```
Original sequence A*34:01:02: *** *** *** ... GCT CCC ...
Replaced sequence A*34:01:02: ATG GCC *TC ... GCT CCC ...
```

#### Generate LILAC sequences

These files are converted to the LILAC allele sequences files with the following command:
```
java -cp lilac.jar com.hartwig.hmftools.lilac.utils.GenerateReferenceSequences \
   -resource_dir /dir_containing/*.txt \ 
   -output_dir /output_dir/
```

## Algorithm

The starting data for the LILAC algorithm is:
- HLA-A, HLA-B, HLA-C and HLA-Y allele sequences from the IMGT/HLA database
- All fragments aligned to HLA-A, HLA-B and HLA-C gene regions

The algorithm has 2 main phases to determine the germline alleles:
1) Elimination phase: Aims to remove allele candidates that are clearly not present
2) Evidence phase: Consider all possible sets of 6 alleles amongst the remaining candidates and chooses the solution that best explains the
   fragments observed

After the germline alleles are determined, LILAC determines the tumor copy number and any somatic mutations in each allele. Note that if
more than 300 bases of the HLA-A,HLA-B and HLA-C coding regions have less than 10 coverage, then LILAC will fail with errors and will not
try to fit the sample.

### Allele set initialisation

The complete set of 4-digit alleles are determined from the 4, 6 and 8 digit alleles in the IMGT/HLA database.

If a 4-digit allele is not explicitly present (e.g. `A*01:02`) but its corresponding 6-digit or 8 digit alleles are
(e.g. `A*01:02:01`, `A*01:02:02`), the numerically lowest is chosen as the representative allele (i.e. `A*01:02:01`).

Some alleles are excluded due to their high similarity closely related pseudogenes, including:
- HLA-H: `A*31:135`, `A*33:191`, `A*02:783`, `B*07:282`
- HLA-Y: `A*30:205`, `A*30:207`, `A*30:225`, `A*30:228`, `A*01:81`, `A*01:237`

### Read pre-processing

All fragments are collected which:
- Are not duplicates
- Have at least 1 read with an alignment overlapping a coding base of HLA-A, HLA-B or HLA-C
- Have all alignments within 1000 bases of an HLA coding region
- Have a mapping quality of at least 1

All reads are trimmed of any bases that overlap the aligned 5' end of its mate (to remove adapter).

### Elimination phase
The elimination phase is primarily an optimization. The goal is simply to reduce the number of possible alleles from ~22k present in the IMGT/HLA database to a manageable number such that the evidence phase can run efficiently.   The principle in the elimination phase is to remove any allele that does not have at least a certain minimal coverage of each of it’s amino acids and bases.    To mitigate the chance of inadvertently eliminating a true allele, common alleles may be recovered at the end of the elimination phase if they have sufficient unique support, but are then penalised in the subsequent evidence phase relative to other candidate alleles.

The steps in the elimination phase are:

#### 1. Nucleotide matrix
At each coding position, create a matrix of (high quality) nucleotide count and determine all bases which are heterozygous across all 6 alleles. 

We do this 3 separate times for HLA-A, B and C. Each time we consider fragments from any of the alignment records that have similar exon boundaries to the type in question. For instance, fragments from the earlier exons which have identical boundaries across all 3 genes will be used to construct the A, B and C matrices, but fragments from the later exons may only contribute to A and B or perhaps only C.  The counts of supporting fragments are then aggregated at each position to construct the nucleotide matrix. 

Fragments with in-frame indels are only included if the indel matches an existing hla type allowing for realignment. Fragments with out-of-frame indels are always excluded (note for the special case of C\*04:09N, a relatively common allele with out of frame indel, it is explicitly rescued at a later stage

During the elimination phase, nucleotide candidates are filtered to include only those with at least max(1,0.000375 \* FragmentCount) high quality (base qual > min(30,medianBaseQuality)) fragment and at least max(2,0.00075 \* FragmentCount) fragments overall support. Subsequently, base quality is not considered.   Sites with more than 1 nucleotide candidate are deemed heterozygous and sites with only 1 are considered to be homozygous across all 6 alleles.

Any alleles with bases that do not match both the heterozygous and homozygous locations of the nucleotide matrix are eliminated.

#### 2. Amino acid matrix
Similarly to the nucleotide matrix, LILAC also constructs a matrix of amino acid candidates.   Again, amino acid candidates are filtered to those where at least max(1,0.000375 \* FragmentCount) fragments support with high base quality (all 3 nucleotides) and at least max(2,0.00075 \* FragmentCount)  fragments over all.  The codon matrix can include inframe insertions and deletions where these match at least one known allele (base quality is not considered).

Exon boundary ‘enrichment’ is applied for all shared amino acids across all 3 genes (amino acids index < 298). This enriches any fragment with nucleotides on one side of an exon boundary with any homozygous nucleotides from the other side so that an amino acid is able to be constructed. 

Similarly to the nucleotide matrix, any alleles with an amino acid or inframe indel that do not match the amino acid matrix are eliminated.

#### 3. Phased haplotypes
In this step we phase the heterozygous amino acid locations and eliminate any alleles that are not supported by phased locations with sufficient overall shared coverage.   

First, we find phased evidence of each consecutive pair of heterozygous codons and record the haplotypes of all fragments containing both codons. This is performed separately for each of HLA-A, HLA-B and HLA-C to account for differences in exon boundaries. Fragments which overlap amino acid 338 onwards will only be assignable to a subset of the alleles, since the exon boundaries differ after this amino acid. The points are only phased if the total coverage is at least 7 fragments per allele [minFragmentsPerAllele] included in the subset at that location (ie. between 14 and 42 fragments with shared coverage depending on the amino acid location). A phased haplotype with only 1 supporting fragment will be removed if the total fragments supporting the pair is 40 or more [minFragmentsToRemoveSingle] (assumed to be a sequencing error).

We then iteratively choose the phased haplotype with the most support and perform the following routine:
- Find other phased evidence that overlaps it.
- Find the minimum number of codon locations required to uniquely identify each phased evidence, eg, if the left evidence has haplotypes [SP, ST] you only need the last codon but if the right evidence has haplotypes [PD, TD, TS] you would need both codon locations.
- Find evidence of fragments that contain all the required codons. As with the paired evidence, there must be at least at least 7 fragments per allele supporting the pair [minFragmentsPerAllele] and a haplotype with only 1 supporting fragment will be removed if the total fragments supporting the pair is 40 or more [minFragmentsToRemoveSingle].
- Check that the new overlapping evidence is consistent with the existing evidence.
- Merge the new evidence with the existing paired evidence.
- Replace the two pieces of used evidence with the new merged evidence.

Once complete, we can eliminate any alleles that do not match the phased evidence. 

#### 4. Recover common alleles
As a fail safe for phasing, any ‘common alleles’ with more than 0.1% population frequency are recovered. The frequencies of alleles are specified in a resource file and are derived from the Hartwig cohort.   

Additionally, C\*04:09N (the most common HLA allele with a frameshift variant) specifically is also rescued if the out of frame indel 6:31237115 CN>C (hg38: chr6:31269338 CN>C) is present.

#### 5. Detect HLA-Y presence 
HLA-Y is a pseudogene that is highly similar to HLA-A and is not present in the human ref genome but is found in approximately 17% of the Hartwig cohort.    The presence of HLA-Y can cause confusion in typing particularly in determining the HLA-A types.    To detect HLA-Y,  LILAC counts the number of fragments that can be assigned uniquely to one of the 3 known HLA-Y alleles and no other candidate alleles.  If at least 1% of fragments align uniquely to HLA-Y then HLA-Y is considered to be present in the sample.   If HLA-Y is found to be present ANY fragment which matches exactly to a HLA-Y allele (uniquely or shared with other alleles) are excluded from further analysis to prevent confusion with highly similar HLA-A alleles. 

#### 6. Test for 2 digit types with unique evidence
To further reduce the number of candidate alleles, If any 2-digit types are sufficiently unique (i.e. uniquely supported by at least 2% of fragments, they are required to contain at least one 4 digit type belonging to that 2 digit type in the evidence phase.   If two 2 digit types from the same gene are found to be sufficiently unique all other alleles are discarded at this point.   If more than two groups are found to be unique the 2 with the highest evidence are supported Any recovered alleles are also discarded at this point unless the 2 digit group has at least one fragment of unique support.   

#### 7. Remove incomplete alleles with insufficient unique evidence
Many alleles in the IMGT database are incomplete (ie contain ‘\*’ characters), all of which are rare in population frequency.  To prevent spurious matches to these wildcard containing alleles in the evidence phase, we eliminate unlikely candidates.  Wildcard containing alleles are eliminated unless they contain at least 2 fragments support for the non wildcard sequence which do not support any remaining candidate allele with a complete sequence defined.

### Evidence phase
In the evidence phase, LILAC evaluates all possible ‘complexes’ (ie. combinations) of remaining alleles that satisfy the following conditions
- There must be either 1 (homozygous) or 2 (heterozygous) alleles belonging to each gene
- At least 1 allele must match each uniquely supported 2 digit type

For each complex, LILAC counts the number of fragments that can be aligned exactly to at least one allele in the complex at all heterozygous locations.  If a fragment has an amino acid which does not match ANY of the amino acid candidates at a heterozygous location, but has at least 1 nucleotide with base qual<min(medianBaseQuality,30), then the amino acid is deemed to match all amino acid candidates   For exon boundaries only exact nucleotide matches are permitted. Any fragments that can be aligned to 2 or more alleles are apportioned equally between the alleles and counted as shared fragments.  

Since many allele definitions have undetermined (‘wildcard’) sequences, particularly in exon 1 and exons 4-8, these require special treatment so that these wildcard alleles are neither unfairly favoured or discriminated against in the fitting.    To achieve this balance, fragments which don’t match an exact sequence in any candidate allele are dropped altogether from consideration such that random sequencing errors or other artefacts overlapping wildcard regions cannot contribute to the complex count for wildcard containing alleles, but any remaining fragments are considered to match an allele if they match all non wildcard sequences (ie. any amino acid is deemed a match to a wildcard).   

Complexes are scored based on the total fragments that can be aligned to at least one allele in the complex, with a small penalty applied base on allele frequency in the population, a penalty for each allele included that was eliminated but subsequently recovered, and a bonus to boost scores of complexes with homozygous allele, and a penalty for solutions with wildcard characters which may cause spurious matches.  The final score is given by:

```
Complex Score = AlignedFragments + FreqPenalty + HomBonus + RecoveryPenalty + Wildcard penalty
where
   FreqPenalty = 0.0018 * SUM[max(log10(Frequency),1e-4)] * AlignedFragments
   HomBonus = 0.0036 * (# of Homozygous alleles) * Fragments
   RecoveryPenalty = 0.0055 * (# of Recovered alleles) * Fragments
   WildcardPenalty = 0.000015 * (# of wildcard characters in alleles) * Fragments
```

If 2 complexes are precisely equally scored, then the solution with the lowest alphabetical 4 digit allele types is chosen. 

The matching unique, apportioned shared, and wildcard (in rare cases where the full allele is not present in the IMGT/HLA database) support for each allele is recorded in both tumor and normal.

2 further performance improvements are made at this step:
- the number of fragments is downsampled to a maximum of 10k for the evaluation.
- if there are predicted to be more than 1 million complexes, then the evidence phase is first performed individually for complexes of 2 alleles per gene to find the top candidates for each of HLA-A, HLA-B & HLA-C.   LILAC retains only the top 5 pairs including each individual allele candidate and then chooses the first 10 unique alleles appearing in the ranked list of pairs, with any common alleles also retained.  The evidence phase is then subsequently run using this reduced set of candidate alleles. 

### Tumor and RNA status of alleles
#### Tumor allele specific copy number
LILAC optionally accepts a tumor bam and a gene copy number file (produced by PURPLE) which contains the minimum copy number and minimum minor allele copy number of each gene.  If a tumor sample is provided, then fragments are counted for each allele in the determined type.   For each of HLA-A, HLA-B & HLA-C, the minor allele copy number is assigned to the allele with the lowest ratio of supporting fragments for the allele in the tumor compared to the normal. The other allele for each is assigned the implied major allele copy number from the gene copy number file.  If a gene is homozygous present in the normal sample, then the minor and major allele copy numbers are arbitrarily assigned in the tumor.

#### Somatic variant assignment to alleles
LILAC optionally accepts a VCF input of somatic small indel and point mutations (called by Sage in the HMF pipeline) and can assign somatic variants to the specific allele which is damaged.

LILAC gathers the set of variants from the vcf (filter = PASS) that overlap either a coding or canonical splice region in any of HLA-A, HLA-B and HLA-C and finds all fragments that contain that variant. LILAC assigns aportions the fragment to the allele which matches the fragment at all heterozygous locations after excluding any somatic variants from the fragment. The allele with the highest matching fragment count is determined to contain the somatic variant. If the variant is assigned to 2 alleles with identical weight it is assigned with 0.5 weight to each. In the case of homozygous alleles, the variant is assigned arbitrarily to the 1st allele.

#### RNA ‘expression’ of alleles
LILAC also optionally accepts a RNA bam.  As per the fragments in the bam are counted for each allele in the determined type.  This can be interpreted as a proxy for allele specific expression.

### QC metrics and PON

LILAC produces a comprehensive set of QC metrics and provides warning statuses if the typing confidence may be diminished. The full list of possible QC warnings is

Warning | Description 
--- | ---
WARN_UNMATCHED_TYPE | all A, B or C types eliminated
WARN_UNMATCHED_INDEL | A novel (non PON) indel is present that was not fit to any allele supported by at least 0.5% of fitted fragments. 
WARN_UNMATCHED_HAPLOTYPE | A novel (non PON) haplotype is present that was not fit to any allele supported by at least 1% of fitted fragments.
WARN_UNMATCHED_AMINO_ACID | A novel (non PON) amino acid is present that was not fit to any allele supported by at least 1% of fitted fragments. 
WARN_UNMATCHED_SOMATIC_VARIANT | A somatic variant from the input vcf could not be assigned to any allele
WARN_LOW_BASE_QUALITY | Median base quality < 25
WARN_LOW_COVERAGE | More than 50 distinct bases of HLA-A, HLA-B and HLA-C combined have less than 10 coverage.
FAIL_LOW_COVERAGE | More than 200 distinct bases of HLA-A, HLA-B and HLA-C combined have less than 10 coverage.

In general the warnings may indicate one of either a novel germline variant/allele, presence of an unusual pseudogene type or a incorrectly typed sample.

In the case of unmatched haplotypes and indels, we find that there are a number of recurrent artefacts found across many samples in our cohort.  We therefore 'PON' filter the following indels and haploytypes prior to calculating QC metrics for the fit

#### INDEL PON
37_Position | 38_Position | Variant | Curation
--- | --- | --- | ---
29912028 | 29944251 | AN>A | Recurrent HLA-H artefact
29911899 | 29944122 | A>ACC | Recurrent HLA-Y artefact
29911318 | 29943541 | CN>C | Recurrent artefact
31324524 | 31356747 | GNN>G | Phased alignment issue
31324528 | 31356751 | G>GGT | Phased alignment issue
29910716 | 29942939 | CN>C | Recurrent artefact
29910657 | 29942880 | G>GT | Recurrent artefact
29912392 | 29944615 | A>AGG | Recurrent artefact
29911899 | 29944122 | A>AC | Recurrent artefact

#### HAPLOTYPE PON
Start | Haplotype | Curation
--- | --- | ---
1 | AVVAPRTLLLLLSGALALTQTWAG | HLA-Y
25 | SHSMRYFSTSVSRPGSGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWMEQEEPEYWDRQTEISKTNAQIDLESLRIALRYYNQSED | HLA-Y
47 | AVGYVDDTQFVRFDSDAASQRMEPRAPWMEQEEPEYWDRQTQISKTNAQIDLESLRIALR | HLA-Y
96 | IDLESLRIALRYYNQSEA | HLA-Y
117 | TIQRMSGCDVGSDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITQRKWEAARQAEQLRAYLEGECMEWLRRYLENGKETLQRT | HLA-Y
206 | DAPKTHMTHHAVSDNEATLRCXALSFYPAEITLTWQR | HLA-Y
206 | DPPKTHMTHYPISDHEATLRCWALG | HLA-H
206 | DPPKTHMTHHPISDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDG | HLA-H
206 | DPPKTHMTHHPISDHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGDG | Likely HLA-Y
242 | RDGEDQTQDTELVETRPAGDGTFQKWASVVVPSGQEQRYTCHVQHEGLPKPLTLRWE | Likely HLA-Y
262 | GIFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWE | HLA-Y
298 | EPSSQPTIPIVGILAGLVLFGAVIAGAVVAAVMWRRKS | Likely HLA-Y
299 | PSSHPTIPIVGILAGLVLFGAVIAGAVVAAVMWRRKS | HLA-Y
299 | IPNLGIVSGPAVLAVLAVLAVLAV | Alignment issue with C\*17 indel
337 | DRKGGSYSQAAS | HLA-Y or HLA-H
351 | IAQGSDVSLTAC | HLA-Y


## Output

The following files are written:

### Solution Summary

Field | Description 
--- | ---
Allele | Allele ID
RefTotal | Total assigned fragments from reference BAM
RefUnique | Fragments uniquely assigned to allele
RefShared | Fragments assigned to allele and others in this solution
RefWild | Fragments matched to a wildcard allele
TumorCopyNumber | Copy number from Tumor/Ref fragment ratio and Purple copy number
TumorTotal | As above for tumor BAM
TumorUnique | As above for tumor BAM
TumorShared | As above for tumor BAM
TumorWild | As above for tumor BAM
RnaTotal | As above for RNA BAM
RnaUnique | As above for RNA BAM
RnaShared | As above for RNA BAM
RnaWild | As above for RNA BAM
SomaticMissense | Matched missense variants
SomaticNonsenseOrFrameshift | Matched nonsense or frameshift variants
SomaticSplice | Matched splice variants
SomaticSynonymous | Matched synonymous variants
SomaticInframeIndel | Matched inframe indels

### QC Metrics

Field | Description 
--- | ---
Status | Either PASS or 1 or more warnings (see below for warning descriptions)
HlaY | HLA-Y allele detected if present, else ‘NONE’
ScoreMargin | Difference in score to second-top solution
NextSolutionAlleles | Allele difference in second-top solution
MedianBaseQuality | Median base quality across all coding bases from all fragments
DiscardedIndels | Discarded fragments due to unknown INDELs
DiscardedIndelMaxFrags | Maximum fragment support for an indel detected but not present in any known allele. 
DiscardedAlignmentFragments | Fragments discarded because 1 read aligns more than 1000 bases from a HLA A,B or C gene
A_LowCoverageBases | Number of bases with less than 15-depth coverage across all coding bases, also for B and C genes
ATypes | # of distinct HLA-A alleles fitted (0,1 or 2)
BTypes | # of distinct HLA-B alleles fitted (0,1 or 2)
CTypes | # of distinct HLA-C alleles fitted (0,1 or 2)
TotalFragments | Total of fragments overlapping a coding base of HLA-A, HLA-B or HLA-C with MAPQ >=1
FittedFragments | Total fragments assigned to fitted alleles (ie. exact or wildcard match at every heterozygous location)
UnmatchedFragments | Fragment that do not match any allele exactly in solution at all heterozygous location (eg. due to sequencing error or mapping error)
UninformativeFragments | Fragments that do not overlap any heterozygous location considered in evidence phase
HlaYFragments | Fragments excluded from fit due to allocation to assignment to HLA-Y.
PercentUnique | Percentage of fitted fragments that are uniquely assigned to 1 allele
PercentShared | Percentage of fitted fragments allocated across multiple alleles
PercentWildcard | Percentage of fitted fragments uniquely assigned to wildcard regions
UnusedAminoAcids | # of amino acids detected with at least 3 fragments support but not present in final allele solution set (excluding PON filtered and regions overlapping wildcards in selected alleles)
UnusedAminoAcidMaxFrags | Maximum fragments support for an unmatched amino acid not present in the final allele solution set
UnusedHaplotypes | # of haplotypes phased with at least 3 fragments support but not present in final allele solution set (excluding PON filtered and regions overlapping wildcards in selected alleles)
UnusedHaplotypeMaxFrags | Maximum support for an unmatched haplotype not present in final allele solution set
SomaticVariantsMatched | Somatic variants supported by solution allele
SomaticVariantsUnmatched | Somatic variants not supported by solution allele


#### Additional output files

File | Description
--- | ---
SAMPLE_ID.lilac.somatic.vcf.gz | Annotation of supplied somatic variants if somatic VCF provided with assigned allele
SAMPLE_ID.fragments.csv | Read details for all BAM fragments
SAMPLE_ID.candidate.coverage.csv | Coverage for all candidate solutions within X% of the top solution's score 
SAMPLE_ID.candidate.fragments.csv | Allocation of each fragment to one or more solutions and which alleles they support 
SAMPLE_ID.HLA-A.aminoacids.txt|Fragment support for each amino acid by HLA gene
SAMPLE_ID.HLA-A.nucleotides.txt|Fragment support for each nucleotide by HLA gene
SAMPLE_ID.candidates.aminoacids.txt|Fragment support for amino acids in the candidate alleles
SAMPLE_ID.candidates.nucleotides.txt|Fragment support for nucleotides in the candidate alleles

The log file contains additional detailed information about the fit including details of all unmatched haplotypes and indels.

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

## Version History and Download Links
- [1.5](https://github.com/hartwigmedical/hmftools/releases/tag/lilac-v1.5.2)
- [1.4](https://github.com/hartwigmedical/hmftools/releases/tag/lilac-v1.4.2)
- [1.3](https://github.com/hartwigmedical/hmftools/releases/tag/lilac-v1.3)
- [1.2](https://github.com/hartwigmedical/hmftools/releases/tag/lilac-v1.2)
- [1.1](https://github.com/hartwigmedical/hmftools/releases/tag/lilac-v1.1)
- [1.0](https://github.com/hartwigmedical/hmftools/releases/tag/lilac-v1.0)

 
