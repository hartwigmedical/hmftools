Running the Hartwig pipeline with Oncoanalyser
==============================================

# Table of contents

- [Overview](#overview)
- [Supported workflows](#supported-workflows)
- [Usage](#usage)
- [Outputs](#outputs)
- [Future improvements](#future-improvements)
- [Acknowledgements](#acknowledgements)

# Introduction
Oncoanalyser ([GitHub](https://github.com/nf-core/oncoanalyser), [nf-core](https://nf-co.re/oncoanalyser/latest/)) is a [Nextflow](https://www.nextflow.io/) 
implementation of the Hartwig Medical Foundation comprehensive cancer DNA/RNA analysis pipeline. In addition to calling simple/structural 
variants and quantifying RNA transcripts, the pipeline also performs downstream analyses such as detecting driver mutations and gene fusions, 
HLA typing, determining tissue of origin, and more. Details on each tool used in the pipeline can be found at [github.com/hartwigmedical/hmftools](https://github.com/hartwigmedical/hmftools).

The pipeline starts from BAM or FASTQ files, and supports data from below modes:

| Data type | Sequencing method                                                                | Paired tumor/normal | tumor-only         |
|-----------|----------------------------------------------------------------------------------|---------------------|--------------------|
| DNA       | Whole genome sequencing (WGS)                                                    | :white_check_mark:  | :white_check_mark: |
| DNA       | Targeted sequencing:<br/> - Whole exome sequencing (WES)<br/> - Panel sequencing | :white_check_mark:  | :white_check_mark: |
| RNA       | Whole transcriptome sequencing (WTS)                                             | -                   | :white_check_mark: |

# Quick start

## 1. Install Nextflow
See: https://www.nextflow.io/docs/latest/install.html

## 2. Install Docker
See: https://docs.docker.com/engine/install/

## 3. Set up resource files

Links:

| Description                                                       | GRCh37                                                                                                                                                                                       | GRCh38                                                                                                                                                                                                       |
|-------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [hmftools](https://github.com/hartwigmedical/hmftools/) resources | [hmf_pipeline_resources.37_v6.0--2.tar.gz](https://storage.googleapis.com/hmf-public/HMFtools-Resources/oncoanalyser/v6_0/37/hmf_pipeline_resources.37_v6.0--2.tar.gz)                       | [hmf_pipeline_resources.38_v6.0--2.tar.gz](https://storage.googleapis.com/hmf-public/HMFtools-Resources/oncoanalyser/v6_0/38/hmf_pipeline_resources.38_v6.0--2.tar.gz)                                       |
| Reference genome                                                  | [Homo_sapiens.GRCh37.GATK.illumina.fasta](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/Homo_sapiens.GRCh37.GATK.illumina.fasta)                               | [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.0/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna)                               |
| Reference genome index                                            | [Homo_sapiens.GRCh37.GATK.illumina.fasta.fai](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/samtools_index/1.16/Homo_sapiens.GRCh37.GATK.illumina.fasta.fai)   | [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.0/samtools_index/1.16/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai)   |
| Reference genome sequence dictionary                              | [Homo_sapiens.GRCh37.GATK.illumina.fasta.dict](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/samtools_index/1.16/Homo_sapiens.GRCh37.GATK.illumina.fasta.dict) | [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.dict](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.0/samtools_index/1.16/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.dict) |
| [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) index            | [bwa-mem2_index/2.13.2](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.1/bwa-mem2_index/2.2.1.tar.gz)                                                             | [bwa-mem2_index/2.13.2](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.1/bwa-mem2_index/2.2.1.tar.gz)                                                                             |
| [GRIDSS](https://github.com/PapenfussLab/gridss) index            | [gridss_index/2.2.1](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.1/gridss_index/2.13.2.tar.gz)                                                                 | [gridss_index/2.2.1](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.1/gridss_index/2.13.2.tar.gz)                                                                                 |
| STAR index (**only required for RNA data**)                       | [star_index/gencode_19/2.7.3a](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/star_index/gencode_19/2.7.3a.tar.gz)                                              | [star_index/gencode_19/2.7.3a](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.0/star_index/gencode_38/2.7.3a.tar.gz)                                                              |

### 3a. Download and extract files
The below commands will set up the resources for human reference genome GRCh37. Adjust the paths and URLs accordingly for GRCh38.

Create a directory structure for the resources:
```shell
mkdir -p $HOME/oncoanalyser/GRCh37/{genome,hmftools}
```

Get genome:
```shell
cd $HOME/oncoanalyser/GRCh37/genome/

wget https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/Homo_sapiens.GRCh37.GATK.illumina.fasta
wget https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/samtools_index/1.16/Homo_sapiens.GRCh37.GATK.illumina.fasta.fai
wget https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/samtools_index/1.16/Homo_sapiens.GRCh37.GATK.illumina.fasta.dict

wget -O bwa-mem2_index.tar.gz https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.1/bwa-mem2_index/2.2.1.tar.gz
tar xvzf bwa-mem2_index.tar.gz --one-top-level

wget -O gridss_index.tar.gz https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.1/gridss_index/2.13.2.tar.gz
tar xvzf gridss_index.tar.gz --one-top-level
```

Get hmftools resources:
```shell
cd $HOME/oncoanalyser/GRCh37/hmftools/
wget https://storage.googleapis.com/hmf-public/HMFtools-Resources/oncoanalyser/v6_0/37/hmf_pipeline_resources.37_v6.0--2.tar.gz
tar xvzf hmf_pipeline_resources.37_v6.0--2.tar.gz --strip-components=1
```

### 3b. Set up config

Create a file called `oncoanalyser_resources.GRCh37.config` which points to the resource paths: 
```
params {
	genomes {
		'GRCh37_hmf' {
			fasta         = "$HOME/oncoanalyser/GRCh37/genome/Homo_sapiens.GRCh37.GATK.illumina.fasta"
			fai           = "$HOME/oncoanalyser/GRCh37/genome/Homo_sapiens.GRCh37.GATK.illumina.fasta.fai"
			dict          = "$HOME/oncoanalyser/GRCh37/genome/Homo_sapiens.GRCh37.GATK.illumina.fasta.dict"
			bwamem2_index = "$HOME/oncoanalyser/GRCh37/genome/bwa-mem2_index/2.13.2/"
			gridss_index  = "$HOME/oncoanalyser/GRCh37/genome/gridss_index/2.2.1/"
		}
	}
	
	ref_data_hmf_data_path = "$HOME/oncoanalyser/GRCh37/hmftools/"
}
```

## 4. Set up sample sheet
Create a file called `samplesheet.csv` which points to the sample inputs:
```
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
COLO829,COLO829,COLO829T,tumor,dna,bam,/path/to/COLO829T.dna.bam
COLO829,COLO829,COLO829R,normal,dna,bam,/path/to/COLO829R.dna.bam
```

BAMs, FASTQs, as well as intermediate outputs can be provided to the sample sheet. Details in [Usage](#usage).

## 6. Run Oncoanalyser with Nextflow
```shell
nextflow run nf-core/oncoanalyser \
-profile docker \
-config oncoanalyser_resources.GRCh37.config \
--mode wgts \
--genome GRCh37_hmf \
--input samplesheet.csv \
--outdir output/
```

## Usage

### Input data

BAM or FASTQ files are the starting inputs for Oncoanalyser.

> [!NOTE] 
> BAM files are expected to meet the below criteria 

All BAM files should be aligned to the Hartwig-distributed
[GRCh37](https://console.cloud.google.com/storage/browser/hmf-public/HMFtools-Resources/ref_genome/37) or
[GRCh38](https://console.cloud.google.com/storage/browser/hmf-public/HMFtools-Resources/ref_genome/38)
reference genomes.

For DNA BAM files: 
- Aligned with [bwa-mem](https://github.com/lh3/bwa), [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) or [DRAGMAP](https://github.com/Illumina/DRAGMAP), with supplementary alignment soft-clipping enabled (i.e. `-Y` argument)

For RNA BAM files:
- Aligned with [STAR](https://github.com/alexdobin/STAR) with some **essential** [settings](../isofox#a-note-on-alignment-and-multi-mapping)
- Duplicates marked with the [Picard](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard)
- Use Ensembl v74 annotations for GRCh37
- Use Ensembl v105 annotations for GRCh38

### Input sample sheet

The sample sheet is a comma separated table where each row represents an input file along with its associated metadata.

| Column        | Description                                                                                         |
| ------------- |-----------------------------------------------------------------------------------------------------|
| group_id      | Group ID for a set of samples and inputs                                                            |
| subject_id    | Subject/patient ID                                                                                  |
| sample_id     | Sample ID                                                                                           |
| sample_type   | Sample type: `tumor`, `normal`                                                                      |
| sequence_type | Sequence type: `dna`, `rna`                                                                         |
| filetype      | File type: `bam`, `fastq`, `bai`, `bam_redux`, etc                                                  |
| info          | For `fastq` file types, specify _library id_ and _lane_, e.g. `library_id:COLO829_library;lane:001` |
| filepath      | Absolute filepath to input file (can be local filepath, URL, S3 URI)                                |

The identifiers provided in the sample sheet are used to set output file paths:

- `group_id`: top-level output directory for analysis files e.g. `output/COLO829/`
- tumor `sample_id`: output prefix for most filenames e.g. `COLO829T.purple.sv.vcf.gz`
- normal `sample_id`: output prefix for some filenames e.g. `COLO829R.cobalt.ratio.pcf`

#### BAM inputs
Below is an example sample sheet with BAM inputs for the whole genome and transcriptome (WGTS) workflow:

```
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
COLO829,COLO829,COLO829T,tumor,dna,bam,/path/to/COLO829T.dna.bam
COLO829,COLO829,COLO829R,normal,dna,bam,/path/to/COLO829R.dna.bam
COLO829,COLO829,COLO829T_RNA,tumor,rna,bam,/path/to/COLO829T.rna.bam
```

> [!NOTE]
> Values in `sample_id` and `filepath` columns must be unique in any sample sheet

> [!NOTE]
> Input filepaths can be absolute local paths, URLs, or S3 URIs

> [!WARNING]
> BAM indexes are expected to exist alongside the respective input BAM but can also be provided as a separate
> samplesheet entry by using the `bai` filetype

#### FASTQ inputs
Below is an example sample sheet with FASTQ inputs for the WGTS workflow:

```
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,info,filepath
COLO829,COLO829,COLO829T,tumor,dna,fastq,library_id:COLO829T_library;lane:001,/path/to/COLO829T.dna.001_R1.fastq.gz;/path/to/COLO829T.dna.001_R2.fastq.gz
COLO829,COLO829,COLO829T,tumor,dna,fastq,library_id:COLO829T_library;lane:002,/path/to/COLO829T.dna.002_R1.fastq.gz;/path/to/COLO829T.dna.002_R2.fastq.gz
COLO829,COLO829,COLO829T,tumor,dna,fastq,library_id:COLO829T_library;lane:003,/path/to/COLO829T.dna.003_R1.fastq.gz;/path/to/COLO829T.dna.003_R2.fastq.gz
COLO829,COLO829,COLO829T,tumor,dna,fastq,library_id:COLO829T_library;lane:004,/path/to/COLO829T.dna.004_R1.fastq.gz;/path/to/COLO829T.dna.004_R2.fastq.gz
COLO829,COLO829,COLO829R,normal,dna,fastq,library_id:COLO829R_library;lane:001,/path/to/COLO829R.dna.001_R1.fastq.gz;/path/to/COLO829R.dna.001_R2.fastq.gz
COLO829,COLO829,COLO829T_RNA,tumor,rna,fastq,library_id:COLO829T_RNA_library;lane:001,/path/to/COLO829T.rna.001_R1.fastq.gz;/path/to/COLO829T.rna.001_R2.fastq.gz
```

The additional `info` column provides the required lane and library info for FASTQ entries with each field delimited by a semicolon.

The forward and reverse FASTQ files are set in the `filepath` column and are also separated by a semicolon, and are _strictly_ ordered 
with forward reads in position one and reverse in position two.

When starting from FASTQ files, reads will be aligned against the selected reference genome using bwa-mem2 (DNA reads) or STAR (RNA reads).

> [!NOTE]
> Only gzipped compressed, non-interleaved pair-end FASTQs are currently supported

#### Run modes

The above examples have provided inputs for the WGTS workflow using paired tumor/normal. However, the below example sample sheets show how
different workflow and/or sample modes can be from BAM files (but also applies to other `sample_type`s e.g. FASTQ files).

Tumor-only DNA:

```
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
COLO829,COLO829,COLO829T,tumor,dna,bam,/path/to/COLO829T.dna.bam
```

Tumor-only DNA and RNA:
```
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
COLO829,COLO829,COLO829T,tumor,dna,bam,/path/to/COLO829T.dna.bam
COLO829,COLO829,COLO829T_RNA,tumor,rna,bam,/path/to/COLO829T.rna.bam
```

RNA only:
```
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
COLO829,COLO829,COLO829T_RNA,tumor,rna,bam,/path/to/COLO829T.rna.bam
```

#### Multiple sample groups
Multiple sample groups can also be provided in a single sample sheet. All rows with the same `group_id` value will be grouped together for 
processing.

```
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
COLO829,COLO829,COLO829T,tumor,dna,bam,/path/to/COLO829T.dna.bam
COLO829,COLO829,COLO829R,normal,dna,bam,/path/to/COLO829R.dna.bam
COLO829,COLO829,COLO829T_RNA,tumor,rna,bam,/path/to/COLO829T.rna.bam
SEQC,SEQC,SEQCT,tumor,dna,bam,/path/to/SEQCT.dna.bam
```

Here the `SEQC` has been added. Since only a tumor DNA BAM is provided for this additional group, just a tumor-only WGS analysis is run 
for the SEQC sample.

#### Inputs other than BAM or FASTQ

It is possible to run Oncoanalyser from any tool or stage as shown in the schematic in [Supported workflows](supported-workflows).

For example, you may already have the inputs data from the WiGiTS pipeline to run [CUPPA](https://github.com/hartwigmedical/hmftools/tree/master/cuppa). 
You would then provide a sample sheet by providing rows with `purple_dir`, `linx_anno_dir` and `isofox_dir` for column `filetype`:

```
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
COLO829,COLO829,COLO829T,tumor,dna,purple_dir,/path/to/purple/dir/
COLO829,COLO829,COLO829T,tumor,dna,linx_anno_dir,/path/to/linx/dir/
COLO829,COLO829,COLO829T,tumor,rna,isofox_dir,/path/to/isofox/dir/
```

Below are all possible values for `filetype`:

- Raw inputs: 
  - `bam`, `bai`,`bam_redux`
  - `fastq`
- Intermediate outputs: 
  - `amber_dir`
  - `bamtools`, `bamtools_dir` 
  - `cobalt_dir`
  - `esvee_vcf`, `esvee_vcf_tbi` 
  - `isofox_dir`
  - `lilac_dir` 
  - `linx_anno_dir` 
  - `pave_vcf`
  - `purple_dir`, 
  - `sage_vcf`, `sage_vcf_tbi`, `sage_append_vcf`
  - `virusinterpreter_dir`
- Running [ORANGE](https://github.com/hartwigmedical/hmftools/tree/master/orange): 
  - `chord_dir` 
  - `sigs_dir`
  - `cuppa_dir`
  - `linx_plot_dir`
  - `sage_dir`

### Example command

To launch oncoanalyser you must provide at least the input samplesheet, the reference genome used for read alignment,
and the desired workflow. When running the targeted sequencing workflow the applicable panel name is also required.

> [!NOTE]
> Setting `-revision` to use a specific version of oncoanalyser is strongly recommended to improve reproducibility and
> stability.

> [!WARNING]
> It is recommended to only run oncoanalyser with Docker, which is done by with `-profile docker`.

#### WGTS workflow command

```
nextflow run nf-core/oncoanalyser \
  -profile docker \
  -revision 0.4.5 \
  --mode wgts \
  --genome GRCh38_hmf \
  --input samplesheet.csv \
  --outdir output/
```

#### Targeted sequencing workflow command

```
nextflow run nf-core/oncoanalyser \
  -profile docker \
  -revision 0.4.5 \
  --mode targeted \
  --panel tso500 \
  --genome GRCh38_hmf \
  --input samplesheet.csv \
  --outdir output/
```

#### Argument descriptions

| Argument       | Group        | Description                                                       |
| -------------- | ------------ | ----------------------------------------------------------------- |
| `-profile`     | Nextflow     | Profile name: `docker` (no other profiles supported at this time) |
| `-revision`    | Nextflow     | Specific oncoanalyser version to run                              |
| `-resume`      | Nextflow     | Use cache from existing run to resume                             |
| `--input`      | oncoanalyser | Samplesheet filepath                                              |
| `--outdir`     | oncoanalyser | Output directory path                                             |
| `--mode`       | oncoanalyser | Workflow name: `wgts`, `targeted`                                 |
| `--panel`      | oncoanalyser | Panel name (only applicable with `--mode targeted`): `tso500`     |
| `--genome`     | oncoanalyser | Reference genome: `GRCh37_hmf`, `GRCh38_hmf`                      |
| `--max_cpus`   | oncoanalyser | Enforce an upper limit of CPUs each process can use               |
| `--max_memory` | oncoanalyser | Enforce an upper limit of memory available to each process        |

## Outputs

The selected results files are written to the output directory and arranged into their corresponding groups by
directories named with the respective `group_id` value from the input samplesheet. Within each group directory, outputs
are further organised by tool.

All intermediate files used by each process are kept in the Nextflow work directory (default: `work/`). Once an analysis
has completed this directory can be removed.

### Sample reports

| Report | Path                                             | Description                                           |
| ------ | ------------------------------------------------ | ----------------------------------------------------- |
| ORANGE | `<group_id>/orange/<tumor_sample_id>.orange.pdf` | PDF summary report of key finding of the HMF pipeline |
| LINX   | `<group_id>/linx/MDX210176_linx.html`            | Interactive HMTL report of all SV plots               |

### Pipeline reports

| Report    | Path                                      | Description                                                        |
| --------- | ----------------------------------------- | ------------------------------------------------------------------ |
| Execution | `pipeline_info/execution_report_*.html`   | HTML report of execution metrics and details                       |
| Timeline  | `pipeline_info/execution_timeline_*.html` | Timeline diagram showing process execution (start/duration/finish) |

## Future Improvements

The following improvements are planned for the next few releases:

- longitudinal analysis of patient samples including ctDNA samples
- cloud-specific instructions and optimisations (ie for AWS, Azure and GCP)

## Acknowledgements

The oncoanalyser pipeline was written by Stephen Watts at the [University of Melbourne Centre for Cancer
Research](https://mdhs.unimelb.edu.au/centre-for-cancer-research) with the support of Oliver Hofmann and the Hartwig
Medical Foundation Australia.
