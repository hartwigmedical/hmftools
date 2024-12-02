<h1>
  <img src="nf-core-oncoanalyser_logo_light.png">
</h1>

Oncoanalyser (links: [GitHub](https://github.com/nf-core/oncoanalyser), [nf-core](https://nf-co.re/oncoanalyser/latest/)) is a 
[Nextflow](https://www.nextflow.io/) a comprehensive cancer DNA/RNA analysis pipeline which runs the tools from the Hartwig Medical 
Foundation ([HMFtools](https://github.com/hartwigmedical/hmftools)).

Please jump to section **[Quick start](#quick-start)** to get started.

### Inputs
Oncoanalyser starts from BAM or FASTQ files, and supports the following modes:

| Data type | Sequencing method                                                                | Paired tumor/normal | Tumor-only         |
|-----------|----------------------------------------------------------------------------------|---------------------|--------------------|
| DNA       | Whole genome sequencing (WGS)                                                    | :white_check_mark:  | :white_check_mark: |
| DNA       | Targeted sequencing:<br/> - Whole exome sequencing (WES)<br/> - Panel sequencing | :white_check_mark:  | :white_check_mark: |
| RNA       | Whole transcriptome sequencing (WTS)                                             | -                   | :white_check_mark: |

### Components
At its core the pipeline performs simple and structural variant calling, but also performs downstream analyses such as detecting driver 
mutations, HLA typing, determining tissue of origin, and more. The table below summarizes the HMFtools Oncoanalyser currently
supports (tools are listed in (roughly) the order they are run:

| Tool                                                                                                                                                       | Description                                                     |
|------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------|
| [Redux](https://github.com/hartwigmedical/hmftools/tree/master/redux)                                                                                      | Mark duplicate reads                                            |
| [Sage](https://github.com/hartwigmedical/hmftools/tree/master/sage)                                                                                        | SNV/indel calling                                               |
| [Pave](https://github.com/hartwigmedical/hmftools/tree/master/pave)                                                                                        | SNV/indel interpretation                                        |
| [Esvee](https://github.com/hartwigmedical/hmftools/tree/master/esvee)                                                                                      | Structural variant calling                                      |
| [Cobalt](https://github.com/hartwigmedical/hmftools/tree/master/cobalt)                                                                                    | Calculate read depth ratios                                     |
| [Amber](https://github.com/hartwigmedical/hmftools/tree/master/amber)                                                                                      | Calculate B-allele frequencies                                  |
| [Purple](https://github.com/hartwigmedical/hmftools/tree/master/purple)                                                                                    | Estimate copy number, purity and ploidy. Identify driver events |
| [Linx](https://github.com/hartwigmedical/hmftools/tree/master/linx)                                                                                        | Structural variant, fusion, and driver event interpretation     |
| [Sigs](https://github.com/hartwigmedical/hmftools/tree/master/sigs)                                                                                        | Quantify mutational signature presence                          |
| [Chord](https://github.com/hartwigmedical/hmftools/tree/master/chord)                                                                                      | Homologous recombination deficiency prediction                  |
| [Cuppa](https://github.com/hartwigmedical/hmftools/tree/master/cuppa)                                                                                      | Tissue of origin prediction                                     |
| [Lilac](https://github.com/hartwigmedical/hmftools/tree/master/lilac)                                                                                      | HLA typing                                                      |
| [Virusbreakend](https://github.com/PapenfussLab/gridss) &<br/>[Virusinterpreter](https://github.com/hartwigmedical/hmftools/tree/master/virus-interpreter) | Viral integration detection	                                    |
| [Isofox](https://github.com/hartwigmedical/hmftools/tree/master/isofox)                                                                                    | Quantify RNA transcripts                                        |
| [Orange](https://github.com/hartwigmedical/hmftools/tree/master/orange)                                                                                    | PDF summary report of all WGS output                            |

## Contents

<!-- TOC -->
    * [Inputs](#inputs)
    * [Components](#components)
  * [Contents](#contents)
  * [Quick start](#quick-start)
    * [1. Install Nextflow](#1-install-nextflow)
    * [2. Install Docker](#2-install-docker)
    * [3. Set up resource files](#3-set-up-resource-files)
      * [3a. Download and extract files](#3a-download-and-extract-files)
      * [3b. Set up config](#3b-set-up-config)
    * [4. Set up sample sheet](#4-set-up-sample-sheet)
    * [5. Run Oncoanalyser with Nextflow](#5-run-oncoanalyser-with-nextflow)
  * [Resource files](#resource-files)
    * [GRCh37](#grch37)
    * [GRCh38](#grch38)
  * [Sample sheet](#sample-sheet)
    * [Sample sheet basic usage](#sample-sheet-basic-usage)
      * [BAM inputs](#bam-inputs)
      * [FASTQ inputs](#fastq-inputs)
      * [Sample modes](#sample-modes)
    * [Sample sheet advanced usage](#sample-sheet-advanced-usage)
      * [Multiple sample groups](#multiple-sample-groups)
      * [Running from REDUX bam](#running-from-redux-bam)
      * [Running specific tools](#running-specific-tools)
  * [Configuration](#configuration)
    * [Arguments](#arguments)
    * [Setting up panel data](#setting-up-panel-data)
    * [Using custom docker images](#using-custom-docker-images)
    * [Using multiple config files](#using-multiple-config-files)
  * [Outputs](#outputs)
    * [Sample reports](#sample-reports)
    * [Pipeline reports](#pipeline-reports)
  * [Future Improvements](#future-improvements)
  * [Acknowledgements](#acknowledgements)
<!-- TOC -->

## Quick start

### 1. Install Nextflow
See: https://www.nextflow.io/docs/latest/install.html

### 2. Install Docker
See: https://docs.docker.com/engine/install/

### 3. Set up resource files

See section **[Resource files](#resource-files)** for URLs to all resources files.

#### 3a. Download and extract files
The below commands will set up the resources for human reference genome GRCh37. 

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
wget https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/bwa_index_image/0.7.17-r1188/Homo_sapiens.GRCh37.GATK.illumina.fasta.img

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

#### 3b. Set up config

Create a file called `oncoanalyser_resources.GRCh37.config` which points to the resource paths: 
```
params {
	genomes {
		'GRCh37_hmf' {
			fasta         = "$HOME/oncoanalyser/GRCh37/genome/Homo_sapiens.GRCh37.GATK.illumina.fasta"
			fai           = "$HOME/oncoanalyser/GRCh37/genome/Homo_sapiens.GRCh37.GATK.illumina.fasta.fai"
			dict          = "$HOME/oncoanalyser/GRCh37/genome/Homo_sapiens.GRCh37.GATK.illumina.fasta.dict"
			img           = "$HOME/oncoanalyser/GRCh37/genome/Homo_sapiens.GRCh37.GATK.illumina.fasta.img"
			bwamem2_index = "$HOME/oncoanalyser/GRCh37/genome/bwa-mem2_index/2.13.2/"
			gridss_index  = "$HOME/oncoanalyser/GRCh37/genome/gridss_index/2.2.1/"
		}
	}
	
	ref_data_hmf_data_path = "$HOME/oncoanalyser/GRCh37/hmftools/"
}
```

### 4. Set up sample sheet
Create a file called `samplesheet.csv` which points to the sample inputs:
```
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
COLO829,COLO829,COLO829T,tumor,dna,bam,/path/to/COLO829T.dna.bam
COLO829,COLO829,COLO829R,normal,dna,bam,/path/to/COLO829R.dna.bam
```

See section **[Sample sheet](#sample-sheet)** for details.

### 5. Run Oncoanalyser with Nextflow
```shell
nextflow run nf-core/oncoanalyser \
-profile docker \
-revision pipeline_v6.0 \
-config oncoanalyser_resources.GRCh37.config \
--mode wgts \
--genome GRCh37_hmf \
--input samplesheet.csv \
--outdir output/
```

## Resource files

### GRCh37

| Type         | Description          | Name                                                                                                                                                                                                |
|--------------|----------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| HMFTools     | HMFTools resources   | [hmf_pipeline_resources.37_v6.0--2.tar.gz](https://storage.googleapis.com/hmf-public/HMFtools-Resources/oncoanalyser/v6_0/37/hmf_pipeline_resources.37_v6.0--2.tar.gz)                              |
| Genome       | FASTA                | [Homo_sapiens.GRCh37.GATK.illumina.fasta](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/Homo_sapiens.GRCh37.GATK.illumina.fasta)                                      |
| Genome       | FASTA index          | [Homo_sapiens.GRCh37.GATK.illumina.fasta.fai](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/samtools_index/1.16/Homo_sapiens.GRCh37.GATK.illumina.fasta.fai)          |
| Genome       | FASTA seq dictionary | [Homo_sapiens.GRCh37.GATK.illumina.fasta.dict](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/samtools_index/1.16/Homo_sapiens.GRCh37.GATK.illumina.fasta.dict)        |
| Genome       | bwa-mem2 index image | [Homo_sapiens.GRCh37.GATK.illumina.fasta.img](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/bwa_index_image/0.7.17-r1188/Homo_sapiens.GRCh37.GATK.illumina.fasta.img) |
| Genome       | bwa-mem2 index       | [bwa-mem2_index/2.2.1.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.1/bwa-mem2_index/2.2.1.tar.gz)                                                              |
| Genome       | GRIDSS index         | [gridss_index/2.13.2.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.1/gridss_index/2.13.2.tar.gz)                                                                |
| Genome (RNA) | STAR index           | [star_index/gencode_19/2.7.3a.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/star_index/gencode_19/2.7.3a.tar.gz)                                              |
| Panel        | TSO500 data          | [panels/tso500_5.34_37--1.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/panels/tso500_5.34_37--1.tar.gz)                                                           |


### GRCh38

| Type         | Description          | Name                                                                                                                                                                                                                |
|--------------|----------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| HMFTools     | HMFTools resources   | [hmf_pipeline_resources.38_v6.0--2.tar.gz](https://storage.googleapis.com/hmf-public/HMFtools-Resources/oncoanalyser/v6_0/38/hmf_pipeline_resources.38_v6.0--2.tar.gz)                                              |
| Genome       | FASTA                | [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.0/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna)                                      |
| Genome       | FASTA index          | [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.0/samtools_index/1.16/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai)          |
| Genome       | FASTA seq dictionary | [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.dict](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.0/samtools_index/1.16/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.dict)        |
| Genome       | bwa-mem2 index image | [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.img](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.0/bwa_index_image/0.7.17-r1188/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.img) |
| Genome       | bwa-mem2 index       | [bwa-mem2_index/2.2.1.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.1/bwa-mem2_index/2.2.1.tar.gz)                                                                              |
| Genome       | GRIDSS index         | [gridss_index/2.13.2.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.1/gridss_index/2.13.2.tar.gz)                                                                                |
| Genome (RNA) | STAR index           | [star_index/gencode_38/2.7.3a.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.0/star_index/gencode_38/2.7.3a.tar.gz)                                                              |
| Panel        | TSO500 data          | [panels/tso500_5.34_38--1.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/panels/tso500_5.34_38--1.tar.gz)                                                                           |

## Sample sheet

The sample sheet is a comma separated table where each row represents an input file along with its associated metadata.

| Column        | Description                                                                                                     |
|---------------|-----------------------------------------------------------------------------------------------------------------|
| group_id      | Group ID for a set of samples and inputs                                                                        |
| subject_id    | Subject/patient ID                                                                                              |
| sample_id     | Sample ID                                                                                                       |
| sample_type   | Sample type: `tumor`, `normal`                                                                                  |
| sequence_type | Sequence type: `dna`, `rna`                                                                                     |
| filetype      | File type: `bam`, `fastq`, etc. See **[Running specific tools](#running-specific-tools)** for all valid values. |
| info          | For `fastq` file types, specify library_id and lane, e.g. `library_id:COLO829_library;lane:001`                 |
| filepath      | Absolute filepath to input file. Can be local filepath, URL, or S3 URI                                          |

### Sample sheet basic usage

#### BAM inputs
Below is an example sample sheet with BAM inputs for the whole genome and transcriptome (WGTS) workflow:

```
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
COLO829,COLO829,COLO829T,tumor,dna,bam,/path/to/COLO829T.dna.bam
COLO829,COLO829,COLO829R,normal,dna,bam,/path/to/COLO829R.dna.bam
COLO829,COLO829,COLO829T_RNA,tumor,rna,bam,/path/to/COLO829T.rna.bam
```

> [!NOTE]
> BAM indexes (.bai files) are expected to exist alongside the respective input BAM unless provided as a separate sample sheet entry by using the `bai` filetype


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

#### Sample modes

The above examples have provided inputs for the WGTS workflow using paired tumor/normal samples. However, the below example sample sheets 
show how different workflow and/or sample modes can be from BAM files (but also applies to other `sample_type`s e.g. FASTQ files).

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

### Sample sheet advanced usage

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

Here the `SEQC` sample has been added. Since only a tumor DNA BAM is provided for this additional group, just a tumor-only WGS analysis is
run for the SEQC sample.

> [!NOTE]
> It is recommended to use one sample sheet per sample group so that errors can easily be isolated.

#### Running from REDUX bam
_TODO_

#### Running specific tools
It is possible to run Oncoanalyser from any [supported tool](#components). For example, you may want to run 
[CUPPA](https://github.com/hartwigmedical/hmftools/tree/master/cuppa) and already have the outputs from [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purple), 
[LINX](https://github.com/hartwigmedical/hmftools/tree/master/linx), and [Virusinterpreter](https://github.com/hartwigmedical/hmftools/tree/master/virus-interpreter). 
In this case, you would provide the outputs from those tools to the sample sheet, specifying entries where `filetype` is `purple_dir`, 
`linx_anno_dir`, and `virusinterpreter_dir`:

```
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
COLO829,COLO829,COLO829T,tumor,dna,purple_dir,/path/to/purple/dir/
COLO829,COLO829,COLO829T,tumor,dna,linx_anno_dir,/path/to/linx/dir/
COLO829,COLO829,COLO829T,tumor,dna,virusinterpreter_dir,/path/to/virus/dir/
```

Please see the respective tool [readmes](#components) for details on which input data is required. 

Below are all valid values for `filetype`:
- Raw inputs: `bam`, `bai`, `fastq`
- REDUX output: `bam_redux`, `redux_jitter_tsv`, `redux_ms_tsv`
- Other tool outputs:
  - `amber_dir`
  - `bamtools`
  - `bamtools_dir`
  - `cobalt_dir`
  - `esvee_vcf`
  - `esvee_vcf_tbi`
  - `isofox_dir`
  - `lilac_dir`
  - `linx_anno_dir`
  - `pave_vcf`
  - `purple_dir`
  - `sage_vcf`
  - `sage_vcf_tbi`
  - `sage_append_vcf`
  - `virusinterpreter_dir`
- ORANGE inputs: `chord_dir`, `sigs_dir`, `cuppa_dir`, `linx_plot_dir`, `sage_dir`

## Configuration

### Arguments

| Argument       | Group        | Description                                                       |
|----------------|--------------|-------------------------------------------------------------------|
| `-profile`     | Nextflow     | Profile name: `docker` (no other profiles supported at this time) |
| `-revision`    | Nextflow     | Specific oncoanalyser version to run                              |
| `-config`      | Nextflow     | Configuration file                                                |
| `-resume`      | Nextflow     | Use cache from existing run to resume                             |
| `--input`      | oncoanalyser | Samplesheet filepath                                              |
| `--outdir`     | oncoanalyser | Output directory path                                             |
| `--mode`       | oncoanalyser | Workflow name: `wgts`, `targeted`                                 |
| `--panel`      | oncoanalyser | Panel name (only applicable with `--mode targeted`): `tso500`     |
| `--genome`     | oncoanalyser | Reference genome: `GRCh37_hmf`, `GRCh38_hmf`                      |
| `--max_cpus`   | oncoanalyser | Enforce an upper limit of CPUs each process can use               |
| `--max_memory` | oncoanalyser | Enforce an upper limit of memory available to each process        |

### Setting up panel data
_TODO_

### Using custom docker images
_TODO_

### Using multiple config files
_TODO_

## Outputs

The selected results files are written to the output directory and arranged into their corresponding groups by
directories named with the respective `group_id` value from the input samplesheet. Within each group directory, outputs
are further organised by tool.

All intermediate files used by each process are kept in the Nextflow work directory (default: `work/`). Once an analysis
has completed this directory can be removed.

### Sample reports

| Report | Path                                             | Description                                           |
|--------|--------------------------------------------------|-------------------------------------------------------|
| ORANGE | `<group_id>/orange/<tumor_sample_id>.orange.pdf` | PDF summary report of key finding of the HMF pipeline |
| LINX   | `<group_id>/linx/MDX210176_linx.html`            | Interactive HMTL report of all SV plots               |

### Pipeline reports

| Report    | Path                                      | Description                                                        |
|-----------|-------------------------------------------|--------------------------------------------------------------------|
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
