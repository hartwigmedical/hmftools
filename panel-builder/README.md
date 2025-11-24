# PanelBuilder

!! Currently in BETA - exact specifications may change as testing is done.

PanelBuilder is a tool to easily create custom panel designs based on simple input features.
You input the genomic features you are interested in and PanelBuilder creates the best set of probes for you.

## Supported Features

- Genes (coding, promoter, UTR, flanks)
- Copy number backbone
- CDR3
- Sample variants
- Custom regions
- Custom structural variants

## Configuration

### Required Arguments

| Argument              | Type | Description                                                      |
|-----------------------|------|------------------------------------------------------------------|
| ref_genome            | Path | Reference genome FASTA file.                                     |
| bwa_index_image       | Path | BWA-MEM index file for the references genome.                    |
| probe_quality_profile | Path | Probe quality resource file. May be GZIP'd.                      |
| output_dir            | Path | Directory in which to store output. Created if it doesn't exist. |

### Optional Arguments

| Argument           | Type                         | Default                     | Description                                                                                                                        |
|--------------------|------------------------------|-----------------------------|------------------------------------------------------------------------------------------------------------------------------------|
| ensembl_data_dir   | Path                         | (none)                      | Ensembl cache directory.                                                                                                           |
| bwa_lib            | Path                         | Search in current directory | Path to BWA-MEM shared library object.                                                                                             |
| genes              | Path                         | (none)                      | Path to TSV file containing desired gene features. If not specified, gene probes are not produced.                                 |
| cn_backbone        | Flag                         | (none)                      | If specified, include copy number backbone probes in the panel.                                                                    |
| cn_backbone_res_kb | Integer                      | 1000                        | Approximate spacing between copy number backbone probes, in kb.                                                                    |
| amber_sites        | Path                         | (none)                      | Path to heterozygous sites TSV file for copy number backbone. May be GZIP'd.                                                       |
| cdr3               | Flag                         | (none)                      | If specified, include CDR3 regions in the panel.                                                                                   |
| sample             | String                       | (none)                      | ID of sample for sample variant probes. If not specified, sample variant probes are not produced.                                  |
| linx_dir           | Path                         | (none)                      | Path to Linx somatic output for sample variant probes.                                                                             |
| linx_germline_dir  | Path                         | (none)                      | Path to Linx germline output for sample variant probes.                                                                            |
| purple_dir         | Path                         | (none)                      | Path to Purple output for sample variant probes.                                                                                   |
| sample_probes      | Integer                      | 500                         | Maximum number of sample variant probes to produce.                                                                                |
| prefer_small_indel | Flag                         | (none)                      | If specified, when selecting nondriver sample variants, prefer variants with lower `IndelLength`.                                  | 
| custom_regions     | Comma-separated list of path | (none)                      | Path(s) to TSV file containing desired custom regions. If not specified, custom region probes are not produced.                    |
| custom_svs         | Path                         | (none)                      | Path to TSV file containing the desired custom structural variants. If not specified, custom structural variants are not produced. |
| threads            | Integer                      | 1                           | Number of threads to use for some parts of the application which support multithreading.                                           |
| output_id          | String                       | (none)                      | Prefix for output files.                                                                                                           |
| verbose_output     | Flag                         | (none)                      | If specified, output more information which may be useful for investigation or debugging. Increases run time.                      |
| log_level          | String                       | `error`                     | `all`/`trace`/`debug`/`info`/`warn`/`error`/`fatal`/`off`                                                                          |

## Example Usage

```shell
java -jar panel-builder.jar \
  -ref_genome Homo_sapiens.GRCh37.GATK.illumina.fasta \
  -ensembl_data_dir ensembl_data_cache/37 \
  -bwa_index_image Homo_sapiens.GRCh37.GATK.illumina.fasta.img \
  -probe_quality_profile panelbuilder_resources/probe_quality_profile.37.tsv.gz \
  -cdr3 \
  -amber_sites panelbuilder_resources/amber_sites.37.tsv.gz \
  -cn_backbone \
  -genes genes_features.tsv \
  -custom_regions custom_regions.tsv \
  -custom_svs custom_svs.tsv \
  -sample COLO829T \
  -linx_dir COLO829T/linx \
  -linx_germline_dir COLO829T/linx_germline \
  -purple_dir COLO829T/purple \
  -output_dir output
```

## Core Concepts

### Target

"Target" refers to a genome region or sequence which we intend to cover with a probe.
Note that since the probe size is fixed, a probe may cover a larger region than the target region.

### Probe Quality Score

A core feature of PanelBuilder is how it intelligently selects probes which we think will give us good results.
To do so, probes are evaluated with a "quality score" (QS), which quantifies how likely a probe is to hybridise with unintended genome region.
If a probe sequence (or part of the sequence) occurs in many places in the genome, then the probe may hybridise there, reducing coverage of the intended region.
For example, if a probe significantly overlaps an Alu region, it is likely to match to the millions of other Alu regions.
The quality score is roughly reciprocal with the number of effective exact matches in the genome.
For example, QS=1 means only one match, the on-target match. QS=0.1 means ten matches, nine of which will be off-target.

The quality score is calculated by aligning the probe sequence with BWA-MEM. The more high-similarity alignments, the lower the quality score.
Since alignment is slow, we have precomputed quality scores across the whole genome, using 40b windows at 20b spacing.
This data is stored in a "probe quality profile" resource file.
The quality score of a probe region can be quickly calculated by taking the minimum of the scores of the windows which overlap (with some interpolation to handle partial overlap).
We have shown that the quality score of a probe calculated from the profile correlates strongly with the quality score calculated from alignment of the full probe.
For cases where the probe is not entirely based off the reference genome (e.g. variants), the application falls back to computing the quality score with alignment of the full probe sequence.

### Accepted/Rejected

All probes produced by PanelBuilder are evaluated against criteria for:

- Quality score
- GC content

(Note the criteria differ between feature types. Some features have no GC constraints while others have strict GC constraints.)

If a probe passes the criteria, then we deem it to be a good probe. The probe is referred to as "accepted" or "acceptable".
Regions which are covered by accepted probes are "accepted" regions.

Probes which do not pass their criteria are referred to as "rejected".
This means we believe these probes do not meet our requirements and likely won't give us good results.
Regions which are not covered by accepted probes are "rejected" regions.

It is important to note that when you request a particular feature, the application may not produce probes which (completely) cover that feature, because they are rejected.
You should check the output to see which probes/regions were rejected, and determine if that is acceptable for you.

## Feature Details & Probe Generation

### Genes

Probes generated to capture different features of genes. Features can be enabled/disabled individually.

Feature methodology:

| Gene feature               | Probe generation                                                                                                                               | Probe evaluation criteria                               |
|----------------------------|------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------|
| Promoter                   | Cover whole region from transcription start to 500b upstream.                                                                                  | `QS>=0.05`                                              |
| UTR                        | One probe covering centre of noncoding exon.                                                                                                   | `QS>=0.05`                                              |
| Coding                     | Cover whole region of coding exons + 10b of splice region.                                                                                     | `QS>=0.05`                                              |
| Exon flanks                | Small intron (3-5kb): one probe in 1kb region centered between probes.<br>Large intron (>5kb): probe in each 1-5kb region 1kb away from exons. | `QS>=0.5`, `0.4<=GC<=0.5`. Select closest to `GC=0.45`  |
| Upstream/downstream flanks | One probe in the 2kb region 1kb upstream/downstream of the gene.                                                                               | `QS>=0.5`, `0.4<=GC<=0.5`. Select closest to `GC=0.45`. |

The canonical transcript of a gene (according to Ensembl) is always included. Additional transcripts may be requested (see input file format below).
If multiple transcripts are requested, they are merged together into a superset of exons from which probes are generated.

#### Gene Feature Input File

TSV file with these columns:

| Column            | Type                 | Description                                                                                |
|-------------------|----------------------|--------------------------------------------------------------------------------------------|
| GeneName          | String               | Ensembl gene name.                                                                         |
| IncludeCoding     | Boolean              | Produce coding exon probes?                                                                |
| IncludeUTR        | Boolean              | Produce UTR probes?                                                                        |
| IncludeExonFlank  | Boolean              | Produce exon flank probes?                                                                 |
| IncludeUpstream   | Boolean              | Produce upstream flank probes?                                                             |
| IncludeDownstream | Boolean              | Produce downstream flank probes?                                                           |
| IncludePromoter   | Boolean              | Produce promoter region probes?                                                            |
| ExtraTransNames   | Comma separated list | Ensembl names of additional transcripts to cover. Canonical transcript is always included. |

Example:

```text
GeneName	IncludeCoding	IncludeUTR	IncludeExonFlank	IncludeUpstream	IncludeDownstream	IncludePromoter	ExtraTransNames
ABCB1	TRUE	FALSE	FALSE	FALSE	FALSE	FALSE	
CDKN2A	TRUE	FALSE	TRUE	TRUE	TRUE	FALSE	ENST00000579755
```

### Copy Number Backbone

A fixed set of probes aimed at improving purity estimation and copy number variation at a genome-wide scale.
These probes are included if you specify the `cn_backbone` argument.

Methodology for all chromosomes except Y:

1. Divide the chromosome into partitions of size `cn_backbone_res_kb`. Exclude 3Mb either side of the centromere.
2. Enumerate Amber heterozygous sites (from `amber_sites` file) within each partition.
3. Filter sites with `0.3<=GNOMAD_FREQ<=0.7`. This ensures the site is likely heterozygous for any individual.
4. Consider all probes centered on the sites.
5. Select the acceptable probe with GC closest to 0.45.

Probe evaluation criteria:

- `QS>=0.8`
- `0.4<=GC<=0.5`.

The methodology for the Y chromosome is similar, except:

- Rather than considering probes on Amber sites, it considers all probes in the partition. This is because there are no heterozygous sites on the Y chromosome.
- QS criteria: `QS=1`.

The strict QS and GC criteria aim to eliminate bias from amplification and hybridisation efficiency, yielding more accurate CN determination.

### CDR3 Regions

A fixed set of probes aimed at capturing all VDJ rearrangements in IG and TCR regions.
These probes are included if you specify the `cdr3` argument.

Methodology:

1. Consider all V and J IG/TCR genes.
2. For each V gene, place one probe aligned to the anchor point, closest to the J gene.
3. For each J gene, place one probe aligned to the anchor point, closest to the V gene.

Probe evaluation criteria:

- `QS>=0.01`.

The QS threshold is low because many of the V/J regions are similar. We have manually determined that this probe set is acceptable.

The list of IG/TCR genes is an embedded resource curated by Hartwig.

With the probes aligned to the anchor of V and J, when recombination occurs, part of the D region will be captured too (which is normally too small to fit a probe).
This allows capturing of the full V+D+J sequence present in the sample.

### Sample Variants

If provided, Linx and Purple output can be used to generate probes covering variants identified in a sample.

Since there may be many variants in a sample, only a subset of variants is selected.
Variants are selected to fill a maximum number of probes, controlled by the `sample_probes` argument.

1. Variants are filtered based on the criteria in the table below. Only variants which pass the filters are selected for probe generation.
2. The selected variants are prioritised. Drivers are given the highest priority, followed by nondrivers.
3. Generate probes from the selected variants, in order of priority, until the maximum number of accepted probes is reached, or there are no more selected variants. Note that a variant may be selected but its probe rejected.

Methodology per variant category:

| Variant category             | Variant selection criteria                                                                                                                                              | Probe evaluation criteria |
|------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------|
| Somatic fusion driver        | (none)                                                                                                                                                                  | `QS>=0.05`                |
| Somatic amplification driver | Highest JCN variant in cluster only                                                                                                                                     | `QS>=0.05`                |
| Somatic deletion driver      | For breakends that flank minimum copy number region                                                                                                                     | `QS>=0.05`                |
| Somatic disruption driver    | `VAF>=0.05`, `AD>=11`                                                                                                                                                   | `QS>=0.05`                |
| Somatic SNV/INDEL driver     | (none)                                                                                                                                                                  | `QS>=0.05`                |
| Germline SV driver           | (none)                                                                                                                                                                  | `QS>=0.05`                |
| Germline SNV/INDEL driver    | (none)                                                                                                                                                                  | `QS>=0.05`                |
| Somatic SNV/INDEL nondriver  | `AD>=11`, `AF>=0.05`, `RC<=3`, `GermlineStatus=DIPLOID`, `IndelLength<=31`. Prioritise lower `IndelLength` (only if configured), then coding, then clonal, then random. | `QS>=0.1`, `0.3<=GC<=0.6` |

Additionally, for all somatic SV drivers, limit variant selection to five breakends per gene.

For each selected variant, one probe is generated that targets the alternate sequence created by the variant. The details are given below:

| Variant type | Probe layout                                                                                                                                                                                                                           |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| SNV          | The reference sequenced is altered by substituting in the variant base. A probe is centered on the resulting sequence.                                                                                                                 |
| INDEL        | The reference sequence is altered by replacing the reference base(s) with the variant bases. A probe is centered on the resulting sequence.                                                                                            |
| SV           | The reference sequence up to the first breakend, then the insert sequence (if any), then the reference sequence from the second breakend. (Breakend orientation is taken into account.) A probe is centered on the resulting sequence. |

### Custom Regions

Arbitrary regions of the reference genome which are completely covered by probes.
Regions are covered with the whole region tiling algorithm described in a subsequent section.

Probe evaluation criteria:

- `QS>=0.1`

#### Custom Regions Input File

TSV file with these columns:

| Column        | Type    | Description                                           |
|---------------|---------|-------------------------------------------------------|
| Chromosome    | String  | Chromosome name as matching the reference genome.     |
| PositionStart | Integer | 1-indexed inclusive start of the region.              |
| PositionEnd   | Integer | 1-indexed inclusive end of the region.                |
| ExtraInfo     | String  | Arbitrary label which will be included in the output. |

Example:
```text
Chromosome	PositionStart	PositionEnd	ExtraInfo
17	7433101	7469631	custom1
1	30429900	30429950	custom2
```

### Custom Structural Variants

Arbitrary structural variants, for which probes are generated in the same manner as described in "Sample Variants".

Probe evaluation criteria:

- `QS>=0.1`

#### Custom Structural Variants Input File

TSV file with these columns:

| Column           | Type        | Description                                                      |
|------------------|-------------|------------------------------------------------------------------|
| ChromosomeStart  | String      | Chromosome name of breakend 1, as matching the reference genome. |
| PositionStart    | Integer     | 1-indexed position of breakend 1.                                |
| OrientationStart | `1` or `-1` | Orientation of breakend 1.                                       |
| ChromosomeEnd    | String      | Chromosome name of breakend 2, as matching the reference genome. |
| PositionEnd      | Integer     | 1-indexed position of breakend 2.                                |
| OrientationEnd   | `1` or `-1` | Orientation of breakend 2.                                       |
| InsertSequence   | String      | Nucleotide sequence inserted in between the two breakends.       |
| ExtraInfo        | String      | Arbitrary label which will be included in the output.            |

Example:

```text
ChromosomeStart	PositionStart	OrientationStart	ChromosomeEnd	PositionEnd	OrientationEnd	InsertSequence	ExtraInfo
1	30510000	1	1	30510020	-1		DEL
1	30600000	-1	1	30600020	1		DUP
1	30700000	1	1	30700020	1		INV
1	30800000	1	2	30800000	-1		TRANS
1	30900000	1	1	30900100	-1	AGGCTGAC	INDEL
```

### Whole Region Tiling

This section describes the algorithm used when a large region is to be fully covered with probes.
This algorithm is applicable for features referring to "cover whole region" or similar wording.

Algorithm:

1. Split the target region into uncovered subregions which do not overlap existing probes.
2. Split the uncovered subregions into acceptable regions where probes may be placed.
This identifies all subregions where overlap with a probe would cause that probe to be rejected according the probe evaluation criteria, and discards those regions.
3. For each acceptable subregion, generate probes to evenly cover ("tile") the subregion. The tiling follows these rules:
    - The tiling always has the minimum number of probes to cover the subregion.
    - Up to 10b on each edge of the subregion may be uncovered before "adding" another probe. This avoids huge probe overlap if the subregion is slightly larger than a multiple of the probe size. (And the probe sequencing is likely to capture these edges anyway).
    - Probes are approximately evenly spaced.
    - The tiling is approximately centered on the subregion.
    - Probes may overlap somewhat.
    - Probes may extend somewhat past the ends of the subregion, as long as this does not place the probe into a rejected region.
    - If the probes cover more bases than the subregion size, the "extra" bases are allocated equally between probe overlap and extension outside the subregion.
4. Subregions of the target region which are not covered by the resulting set of probes (and are not the 10b of allowable uncovered edge) are marked as rejected.

### Probe Overlap Handling

Probes are generated in this order, with probes generated first having priority over subsequent probes:

1. Genes
2. Custom regions
3. Custom structural variants
4. Copy number backbone
5. CDR3
6. Sample variants

In general, if a probe's target region is already covered by an existing probe, the new probe is not included in the panel.

Exceptions:

- Whole region tiling: The algorithm attempts to avoid overlapping existing probes, although it is not guaranteed.
- Variant probes: If the variant introduces at least 5b of difference to the reference genome, the probe is always included, no matter the amount of overlap with existing probes.

## Output

Main outputs:

| File               | Description                                                   |
|--------------------|---------------------------------------------------------------|
| panel_probes.tsv   | Full information for each probe in the panel.                 |
| panel_probes.fasta | Base sequences of probes in the panel.                        |
| rejections.tsv     | Full information for each uncovered region or rejected probe. |

Informational/visualisation/debugging outputs:

| File                     | Description                                                                                                                                            |
|--------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------|
| panel_probes.bed         | Regions of probes which correspond to exact locations in the reference genome. Excludes variant probes.                                                |
| targeted_regions.bed     | Regions which the probes are targeting (subset of probe regions). Includes variant probes but only the parts which correspond to the reference genome. |
| rejections.bed           | Regions which were rejected. Excludes variant probes.                                                                                                  |                                                                                                
| gene_stats.tsv           | Statistics on probes on a per-gene basis. Only produced if gene features were requested.                                                               |
| sample_variant_info.tsv  | Additional information used in processing on a per-variant basis. Only produced if sample variants probes were requested.                              |
| candidate_regions.bed.gz | All regions evaluated for suitability.                                                                                                                 |
| candidate_probes.tsv.gz  | All probes evaluated for suitability. Only produced if `verbose_output` is specified.                                                                  |

All output files will be prefixed by `output_id` if specified.
