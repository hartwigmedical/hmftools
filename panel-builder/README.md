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

## Configuration

### Required Arguments

| Argument              | Type | Description                                                      |
|-----------------------|------|------------------------------------------------------------------|
| ref_genome            | Path | Reference genome FASTA file.                                     |
| ensembl_data_dir      | Path | Ensembl cache directory.                                         |
| bwa_index_image       | Path | BWA-MEM index file for the references genome.                    |
| probe_quality_profile | Path | Probe quality resource file. May be GZIP'd.                      |
| output_dir            | Path | Directory in which to store output. Created if it doesn't exist. |

### Optional Arguments

| Argument          | Type   | Default                     | Description                                                                                                                          |
|-------------------|--------|-----------------------------|--------------------------------------------------------------------------------------------------------------------------------------|
| bwa_lib           | Path   | Search in current directory | Path to BWA-MEM shared library object.                                                                                               |
| genes             | Path   | (none)                      | Path to TSV file containing desired gene features. If not specified, gene probes are not produced.                                   |  
| amber_sites       | Path   | (none)                      | Path to heterozygous sites TSV file for copy number backbone. May be GZIP'd. If not specified, copy number backbone is not produced. |
| cdr3              | Flag   | (none)                      | If specified, include CDR3 regions in the panel.                                                                                     |
| sample            | String | (none)                      | ID of sample for sample variant probes. If not specified, sample variant probes are not produced.                                    |
| linx_dir          | Path   | (none)                      | Path to Linx somatic output for sample variant probes.                                                                               |
| linx_germline_dir | Path   | (none)                      | Path to Linx germline output for sample variant probes.                                                                              |
| purple_dir        | Path   | (none)                      | Path to Purple output for sample variant probes.                                                                                     |
| custom_regions    | Path   | (none)                      | Path to TSV file containing desired custom regions. If not specified, custom region probes are not produced.                         |
| output_id         | String | (none)                      | Prefix for output files.                                                                                                             |
| verbose_output    | Flag   | (none)                      | If specified, output more information which may be useful for investigation or debugging. Increases run time.                        |
| log_level         | String | `error`                     | `all`/`trace`/`debug`/`info`/`warn`/`error`/`fatal`/`off`                                                                            |

## Example Usage

```shell
java -jar panel-builder.jar \
  -ref_genome Homo_sapiens.GRCh37.GATK.illumina.fasta \
  -ensembl_data_dir ensembl_data_cache/37 \
  -bwa_index_image Homo_sapiens.GRCh37.GATK.illumina.fasta.img \
  -probe_quality_profile panelbuilder_resources/probe_quality_profile.37.tsv.gz \
  -cdr3 \
  -amber_sites panelbuilder_resources/amber_sites.37.tsv.gz \
  -genes genes_features.tsv \
  -custom_regions custom_regions.tsv \
  -sample COLO829T \
  -linx_dir COLO829T/linx \
  -linx_germline_dir COLO829T/linx_germline \
  -purple_dir COLO829T/purple \
  -output_dir output
```

## Core Concepts

TODO

- target
- quality score
- accepted/rejected

## Feature Details & Probe Generation

TODO

### Genes

TODO

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

```
GeneName	IncludeCoding	IncludeUTR	IncludeExonFlank	IncludeUpstream	IncludeDownstream	IncludePromoter	ExtraTransNames
ABCB1	TRUE	FALSE	FALSE	FALSE	FALSE	FALSE	
CDKN2A	TRUE	FALSE	TRUE	TRUE	TRUE	FALSE	ENST00000579755
```

### Copy Number Backbone

TODO

### CDR3 Regions

TODO

### Sample Variants

TODO

### Custom Regions

TODO

#### Custom Regions Input File

TSV file with these columns:

| Column        | Type    | Description                                           |
|---------------|---------|-------------------------------------------------------|
| Chromosome    | String  | Chromosome name as matching the reference genome.     |
| PositionStart | Integer | 1-indexed inclusive start of the region.              |
| PositionEnd   | Integer | 1-indexed inclusive end of the region.                |
| ExtraInfo     | String  | Arbitrary label which will be included in the output. |

Example:
```
Chromosome	PositionStart	PositionEnd	ExtraInfo
17	7433101	7469631	custom1
1	30429900	30429950	custom2
```

## Output

TODO
