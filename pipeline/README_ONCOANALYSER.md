<h1>
  <img src="nf-core-oncoanalyser_logo_light.png">
</h1>

Oncoanalyser (links: [GitHub](https://github.com/nf-core/oncoanalyser), [nf-core](https://nf-co.re/oncoanalyser/latest/)) is a 
[Nextflow](https://www.nextflow.io/) implementation of the Hartwig Medical Foundation DNA and RNA sequencing analysis pipeline.

Except for read alignment, the pipeline uses tools from [HMFtools](https://github.com/hartwigmedical/hmftools/tree/master/):
- Read mapping: [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) (DNA), [STAR](https://github.com/alexdobin/STAR) (RNA)
- Read deduplication and unmapping: [REDUX](https://github.com/hartwigmedical/hmftools/tree/master/redux)
- SNV, MNV and INDEL calling: [SAGE](https://github.com/hartwigmedical/hmftools/tree/master/sage), [PAVE](https://github.com/hartwigmedical/hmftools/tree/master/pave)
- CNV calling: [COBALT](https://github.com/hartwigmedical/hmftools/tree/master/cobalt), [AMBER](https://github.com/hartwigmedical/hmftools/tree/master/amber), [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purple)
- SV calling: [ESVEE](https://github.com/hartwigmedical/hmftools/tree/master/esvee)
- SV and driver event interpretation: [LINX](https://github.com/hartwigmedical/hmftools/tree/master/linx)
- Oncoviral detection: [VIRUSbreakend](https://github.com/PapenfussLab/gridss), [VirusInterpreter](https://github.com/hartwigmedical/hmftools/tree/master/virus-interpreter)
- HLA typing: [LILAC](https://github.com/hartwigmedical/hmftools/tree/master/lilac)
- HRD prediction: [CHORD](https://github.com/hartwigmedical/hmftools/tree/master/chord)
- Tissue of origin prediction: [CUPPA](https://github.com/hartwigmedical/hmftools/tree/master/cuppa)
- Mutational signature fitting: [Sigs](https://github.com/hartwigmedical/hmftools/tree/master/sigs)
- RNA transcript quantification: [ISOFOX](https://github.com/hartwigmedical/hmftools/tree/master/isofox)
- Summary report PDF: [ORANGE](https://github.com/hartwigmedical/hmftools/tree/master/orange)

Oncoanalyser supports the following sequencing and sample setups:

| Data type | Sequencing method                                                                | Paired tumor/normal | Tumor-only         |
|-----------|----------------------------------------------------------------------------------|---------------------|--------------------|
| DNA       | Whole genome sequencing (WGS)                                                    | :white_check_mark:  | :white_check_mark: |
| DNA       | Targeted sequencing:<br/> - Whole exome sequencing (WES)<br/> - Panel sequencing | :white_check_mark:  | :white_check_mark: |
| RNA       | Whole transcriptome sequencing (WTS)                                             | -                   | :white_check_mark: |

## Getting started

This section will assume that the analysis starts from tumor/normal BAMs mapped to GRCh37.

### 1. Install Nextflow
See: **https://www.nextflow.io/docs/latest/install.html**

### 2. Install Docker
See: **https://docs.docker.com/engine/install/**

### 3. Set up resource files

Download and extract the reference genome and HMFTools resources using these **[links](#links)**.

Create a file called `oncoanalyser.config` which points to the resource file paths:

```
params {
   genomes {
      'GRCh37_hmf' {
         fasta         = "/path/to/Homo_sapiens.GRCh37.GATK.illumina.fasta"
         fai           = "/path/to/Homo_sapiens.GRCh37.GATK.illumina.fasta.fai"
         dict          = "/path/to/Homo_sapiens.GRCh37.GATK.illumina.fasta.dict"
         img           = "/path/to/Homo_sapiens.GRCh37.GATK.illumina.fasta.img"
         bwamem2_index = "/path/to/bwa-mem2_index/"
         gridss_index  = "/path/to/gridss_index/"
      }
   }

   ref_data_hmf_data_path = "/path/to/hmf_pipeline_resources/"
}
```

> [!TIP]
> Jump to section: **[Configuration files](#configuration-files)**

### 4. Set up sample sheet
Create a file called `sample_sheet.csv` which points to the sample inputs:

```csv
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
COLO829,COLO829,COLO829T,tumor,dna,bam,/path/to/COLO829T.dna.bam
COLO829,COLO829,COLO829R,normal,dna,bam,/path/to/COLO829R.dna.bam
```

BAM and BAI files for the above COLO829 test sample can be downloaded from [here](test_data/).

> [!TIP]
> Jump to section: **[Sample sheet](#sample-sheet)**

### 5. Run Oncoanalyser with Nextflow

```shell
nextflow run nf-core/oncoanalyser \
-profile docker \
-revision pipeline_v6.0 \
-config oncoanalyser.config \
--mode wgts \
--genome GRCh37_hmf \
--input sample_sheet.csv \
--outdir output/
```

> [!TIP]
> Jump to section: **[Command line interface](#command-line-interface--cli-)**

## Table of contents
<!-- TOC -->
  * [Getting started](#getting-started)
    * [1. Install Nextflow](#1-install-nextflow)
    * [2. Install Docker](#2-install-docker)
    * [3. Set up resource files](#3-set-up-resource-files)
    * [4. Set up sample sheet](#4-set-up-sample-sheet)
    * [5. Run Oncoanalyser with Nextflow](#5-run-oncoanalyser-with-nextflow)
  * [Table of contents](#table-of-contents)
  * [Command line interface (CLI)](#command-line-interface--cli-)
    * [Running Oncoanalyser](#running-oncoanalyser)
    * [Nextflow arguments](#nextflow-arguments)
    * [Oncoanalyser arguments](#oncoanalyser-arguments)
    * [Params file](#params-file)
  * [Sample sheet](#sample-sheet)
    * [BAM inputs](#bam-inputs)
    * [FASTQ inputs](#fastq-inputs)
    * [Sample modes](#sample-modes)
    * [Multiple sample groups](#multiple-sample-groups)
    * [Running from REDUX BAM](#running-from-redux-bam)
    * [Running specific tools](#running-specific-tools)
  * [Configuration files](#configuration-files)
    * [Basic config example](#basic-config-example)
    * [Oncoanalyser arguments as config](#oncoanalyser-arguments-as-config)
    * [Multiple config files](#multiple-config-files)
  * [Resource files](#resource-files)
    * [Links](#links)
    * [Configuring WGTS resource files](#configuring-wgts-resource-files)
    * [Configuring panel / targeted sequecing resource files](#configuring-panel--targeted-sequecing-resource-files)
  * [Configuring processes](#configuring-processes)
    * [Compute resources](#compute-resources)
    * [Docker images](#docker-images)
  * [Outputs](#outputs)
    * [Sample reports](#sample-reports)
    * [Pipeline reports](#pipeline-reports)
  * [Acknowledgements](#acknowledgements)
<!-- TOC -->

## Command line interface (CLI)

### Running Oncoanalyser
We use the `nextflow run` command to run the Oncoanalyser. Below is an example that covers some of the useful arguments.
[Nextflow-specific arguments](https://www.nextflow.io/docs/latest/reference/cli.html) have a single hyphen (`-`) while Oncoanalyser-specific
arguments have two hyphens (`--`).

```shell
nextflow run nf-core/oncoanalyser \
-profile docker \
-revision dev \
-config hmf_pipeline_resources.config \
--mode wgts \
--genome GRCh37_hmf \
--input samplesheet.csv \
--outdir output/ \
--max_cpus 32 \
--max_memory 128.GB \
-resume
```

The above command will automatically pull the Oncoanalyser [git repo](https://github.com/nf-core/oncoanalyser). However, we can point
`nextflow run` to a local Oncoanalyser repo (e.g. one we've manually pulled), which can be useful for debugging.
This will run repo with the currently checked out commit and is incompatible with the `-revision` argument.
```shell
nextflow run /path/to/oncoanalyser_repo \
# other arguments
```

### Nextflow arguments
All arguments for `nextflow run` are documented in the [CLI reference](https://www.nextflow.io/docs/latest/reference/cli.html#run). The
below table list the ones that are mandatory or useful.

| Argument&emsp; | Description                                                                                                                                                                                                                            |
|:---------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `-config`      | Path to a [configuration file](#configuration-files). Can be specified [multiple times](#multiple-config-files)                                                                                                                        |
| `-params-file` | Path to a [params file](https://www.nextflow.io/docs/latest/cli.html#pipeline-parameters) containing Oncoanalyser command line arguments                                                                                               |
| `-profile`     | Compute profile. Only `docker` is currently supported.                                                                                                                                                                                 |
| `-revision`    | A specific Oncoanalyser branch/tag to run. See the Oncoanalyser [GitHub](https://github.com/nf-core/oncoanalyser) for available branches/tags                                                                                          |
| `-resume`      | [Resume](https://www.nextflow.io/docs/latest/cache-and-resume.html#work-directory) from cached results (by default the previous run). Useful if you've cancelled a run with `CTRL+C`, or a run has crashed and you've fixed the issue. |
| `-stub`        | Dry run. Under the hood, Oncoanalyser runs `touch <outputfile>` rather than actually running the tools. Useful for testing if the arguments and configuration files provided are correct.                                              |
| `-work-dir`    | Path to a directory where Nextflow will put temporary files for each step in the pipeline. If this is not specified, Nextflow will create the `work/` directory in the current directory                                               |
| `-help`        | Show all Nextflow command line arguments and their descriptions                                                                                                                                                                        |

### Oncoanalyser arguments

| Argument&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; | Description                                                                                                                                                                                                                                                             |
|:---------------------------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--input`                                          | Path to a [sample sheet](#sample-sheet)                                                                                                                                                                                                                                 |
| `--outdir`                                         | Path to the output directory. While a process/tool is running, files are temporarily stored in the work directory (see: `-work-dir` [argument](#nextflow-arguments)). Only when the process completes are the files copied to the output directory.                     |
| `--mode`                                           | Can be:<br/>- `wgts`: Whole genome sequencing and/or whole transcriptome sequencing analysis<br/>- `targeted`: Targeted sequencing analysis (e.g. for panel or whole exome sequencing)                                                                                  |
| `--genome`                                         | Reference genome version. Can be `GRCh37_hmf` or `GRCh38_hmf`                                                                                                                                                                                                           |
| `--ref_data_hmf_data_path`                         | Path to HMFTools resource files                                                                                                                                                                                                                                         |
| `--ref_data_panel_data_path`                       | Path to panel resource files                                                                                                                                                                                                                                            |
| `--ref_data_hla_slice_bed`                         | Path to HLA slice BED file                                                                                                                                                                                                                                              |
| `--email`                                          | Email address for completion summary                                                                                                                                                                                                                                    |
| `--create_stub_placeholders`                       | Create placeholders for resource files during stub run.                                                                                                                                                                                                                 |
| `--max_cpus`                                       | Enforce an upper limit of CPUs each process can use, e.g. `16`                                                                                                                                                                                                          |
| `--max_memory`                                     | Enforce an upper limit of memory available to each process, e.g. `32.GB`                                                                                                                                                                                                |
| `--max_fastq_records`                              | When positive, will use [fastp](https://github.com/OpenGene/fastp) to split fastq files so that each resultant fastq file has no more than max_fastq_records records. When nonpositive, fastp is not used and the provided fastq files are passed as-is to the aligner. |
| `--processes_exclude`<sup>1</sup>                  | A comma separated list specifying which processes to skip (e.g. `--processes_exclude lilac,virusinterpreter`). Note: Downstream processes depending on the output of an upstream tool will also be skipped.                                                             |
| `--processes_manual`                               | Only run processes provided in `--processes_include`                                                                                                                                                                                                                    |
| `--processes_include`                              | When also specifying `--processes_manual`, a comma separated list specifying which processes to include (e.g. `--processes_include lilac,virusinterpreter`).                                                                                                            |
| `--prepare_reference_only`                         | Only stage reference genome and resource files                                                                                                                                                                                                                          |
| `--isofox_read_length`                             | User defined RNA read length used for ISOFOX                                                                                                                                                                                                                            |
| `--isofox_gc_ratios`                               | User defined User defined ISOFOX expected GC ratios file                                                                                                                                                                                                                |
| `--isofox_counts`                                  | User defined ISOFOX expected counts files (read length dependent)                                                                                                                                                                                                       |
| `--isofox_tpm_norm`                                | User defined ISOFOX TPM normalisation file for panel data                                                                                                                                                                                                               |
| `--isofox_gene_ids`                                | User defined ISOFOX gene list file for panel data.                                                                                                                                                                                                                      |
| `--isofox_functions`                               | Semicolon-separated list of ISOFOX functions to run. Default: `TRANSCRIPT_COUNTS;ALT_SPLICE_JUNCTIONS;FUSIONS;RETAINED_INTRONS`                                                                                                                                         |

Notes:
1. Valid process names are: `alignment`, `amber`, `bamtools`, `chord`, `cobalt`, `cuppa`, `esvee`, `isofox`, `lilac`, `linx`, `neo`,
   `orange`, `pave`, `purple`, `redux`, `sage`, `sigs`, `virusinterpreter`

## Sample sheet

The sample sheet is a comma separated table with the following columns:
- `group_id`: Groups `subject_id` and `sample_id` into the same experiment 
- `subject_id`
- `sample_id`
- `sample_type`: `tumor` or `normal`
- `sequence_type`: `dna` or `rna`
- `filetype`: `bam`, `bai`, `fastq`, or see **[Running specific tools](#running-specific-tools)** for other valid values
- `info`: Sequencing library and lane info for **[FASTQ inputs](#fastq-inputs)**
- `filepath`: Absolute filepath to input file. Can be local filepath, URL, or S3 URI

### BAM inputs
Below is an example sample sheet with BAM files for a tumor/normal WGS run:

```csv
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
COLO829,COLO829,COLO829T,tumor,dna,bam,/path/to/COLO829T.dna.bam
COLO829,COLO829,COLO829R,normal,dna,bam,/path/to/COLO829R.dna.bam
```

BAM indexes (.bai files) are expected to be in the same directory as the BAM files. Alternatively, provide the BAM index path by 
providing `bai` under column `filetype`:

```csv
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
COLO829,COLO829,COLO829T,tumor,dna,bam,/path/to/COLO829T.dna.bam
COLO829,COLO829,COLO829T,tumor,dna,bai,/path/to/COLO829T.dna.bam.bai
COLO829,COLO829,COLO829R,normal,dna,bam,/path/to/COLO829R.dna.bam
COLO829,COLO829,COLO829R,normal,dna,bai,/path/to/COLO829R.dna.bam.bai
```

### FASTQ inputs
Below is an example sample sheet with FASTQ files for a tumor/normal WGS run:

```csv
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,info,filepath
COLO829,COLO829,COLO829T,tumor,dna,fastq,library_id:S1;lane:001,/path/to/COLO829T_S1_L001_R1_001.fastq.gz;/path/to/COLO829T_S1_L001_R2_001.fastq.gz
COLO829,COLO829,COLO829T,tumor,dna,fastq,library_id:S1;lane:002,/path/to/COLO829T_S1_L002_R1_001.fastq.gz;/path/to/COLO829T_S1_L002_R2_001.fastq.gz
COLO829,COLO829,COLO829R,normal,dna,fastq,library_id:S2;lane:001,/path/to/COLO829R_S2_L001_R1_001.fastq.gz;/path/to/COLO829R_S2_L001_R2_001.fastq.gz
COLO829,COLO829,COLO829R,normal,dna,fastq,library_id:S2;lane:002,/path/to/COLO829R_S2_L002_R1_002.fastq.gz;/path/to/COLO829R_S2_L002_R2_001.fastq.gz
```

Comments:
- Under `info`, provide the sequencing library and lane info separated by `;`
- Under `filepath`, provide the forward ('R1') and reverse ('R2') FASTQ files separated by `;`

> [!NOTE]
> Only gzip compressed, non-interleaved pair-end FASTQ files are currently supported

### Sample modes

Providing `sample_type` and `sequence_type` in different combinations allows Oncoanalyser to run in different sample modes. The below sample
sheets use BAM files, but different sample modes can also be specified for FASTQ files.

**Tumor-only DNA**

```csv
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
COLO829,COLO829,COLO829T,tumor,dna,bam,/path/to/COLO829T.dna.bam
```

**Tumor-only RNA**
```csv
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
COLO829,COLO829,COLO829T_RNA,tumor,rna,bam,/path/to/COLO829T.rna.bam
```

**Tumor/normal DNA, tumor-only RNA**

```csv
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
COLO829,COLO829,COLO829T,tumor,dna,bam,/path/to/COLO829T.dna.bam
COLO829,COLO829,COLO829R,normal,dna,bam,/path/to/COLO829R.dna.bam
COLO829,COLO829,COLO829_RNA,tumor,dna,bam,/path/to/COLO829R.rna.bam
```

### Multiple sample groups

Multiple sample groups can also be provided in a single sample sheet. All rows with the same `group_id` value will be grouped together for
processing.

```
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
COLO829,COLO829,COLO829T,tumor,dna,bam,/path/to/COLO829T.dna.bam
COLO829,COLO829,COLO829R,normal,dna,bam,/path/to/COLO829R.dna.bam
PD10010,PD10010,PD10010T,tumor,dna,bam,/path/to/PD10010T.dna.bam
PD10010,PD10010,PD10010R,normal,dna,bam,/path/to/PD10010R.dna.bam
```

> [!NOTE]
> It is still recommended to use one sample sheet per sample group so that errors can easily be isolated.

### Running from REDUX BAM
Read mapping with [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) followed by read pre-processing with [REDUX](https://github.com/hartwigmedical/hmftools/tree/master/redux) 
are the pipeline steps that take the most time and compute resources. Thus, we can re-run Oncoanalyser from a REDUX BAM if it is already 
exists, e.g. due to updates to downstream [HMFtools](https://github.com/hartwigmedical/hmftools/tree/master/).

Simply provide the REDUX BAM path, specifying `bam_redux` under `filetype`:

```csv
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
COLO829,COLO829,COLO829T,tumor,dna,bam_redux,/path/to/COLO829T.dna.redux.bam
```

The `*.jitter_params.tsv` and `*.ms_table.tsv.gz` REDUX output files are expected to be in the same directory as the REDUX BAM. If these 
files are located elsewhere, their paths can also be explicitly provided by specifying `redux_jitter_tsv` and `redux_ms_tsv` under `filetype`:

```csv
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
COLO829,COLO829,COLO829T,tumor,dna,bam_redux,/path/to/COLO829T.dna.redux.bam
COLO829,COLO829,COLO829T,tumor,dna,redux_jitter_tsv,/path/to/COLO829T.dna.jitter_params.tsv
COLO829,COLO829,COLO829T,tumor,dna,redux_ms_tsv,/path/to/COLO829T.dna.ms_table.tsv.gz
```

### Running specific tools
It is possible to run Oncoanalyser from any tool from [HMFtools](https://github.com/hartwigmedical/hmftools/tree/master/). For example, you may want to 
run [CUPPA](https://github.com/hartwigmedical/hmftools/tree/master/cuppa) and already have the outputs from
[PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purple), 
[LINX](https://github.com/hartwigmedical/hmftools/tree/master/linx), 
and [VirusInterpreter](https://github.com/hartwigmedical/hmftools/tree/master/virus-interpreter). In this case, you would provide the 
outputs from those tools to the sample sheet, specifying entries where `filetype` is `purple_dir`, `linx_anno_dir`, and 
`virusinterpreter_dir`:

```csv
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
COLO829,COLO829,COLO829T,tumor,dna,purple_dir,/path/to/purple/dir/
COLO829,COLO829,COLO829T,tumor,dna,linx_anno_dir,/path/to/linx/dir/
COLO829,COLO829,COLO829T,tumor,dna,virusinterpreter_dir,/path/to/virus/dir/
```

Please see the respective tool [READMEs](https://github.com/hartwigmedical/hmftools/tree/master/) for details on which input data is required. 

Below are all valid values for `filetype`:

| Type               | Values                                                                                                                                                                                                                           |
|--------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Raw inputs         | `bam`, `bai`, `fastq`                                                                                                                                                                                                            |
| REDUX output       | `bam_redux`, `redux_jitter_tsv`, `redux_ms_tsv`                                                                                                                                                                                  |
| Other tool outputs | `amber_dir`, `bamtools`, `bamtools_dir`, `cobalt_dir`, `esvee_vcf`, `esvee_vcf_tbi`, `isofox_dir`, `lilac_dir`, `linx_anno_dir`, `pave_vcf`, `purple_dir`, `sage_vcf`, `sage_vcf_tbi`, `sage_append_vcf`, `virusinterpreter_dir` |
| ORANGE inputs      | `chord_dir`, `sigs_dir`, `cuppa_dir`, `linx_plot_dir`, `sage_dir`                                                                                                                                                                |

## Configuration files
[Nextflow configuration files](https://www.nextflow.io/docs/latest/config.html) can be used to configure Oncoanalyser. This section 
summarizes concepts of Nextflow configuration files that are relevant for using Oncoanalyser.

For details on specific configurations, please jump to the relevant section:
- [Configuring resource files](#resources-files)
- [Configuring processes](#configuring-processes)

### Basic config example
Config items can be declared using [blocks](#https://www.nextflow.io/docs/latest/config.html#blocks), where curly brackets define the scope 
of the encapsulated config items. The below example has the `params` and `process` scopes, with `workDir` being 
[un-scoped](https://www.nextflow.io/docs/latest/reference/config.html#unscoped-options):

```
workDir = '/path/to/work/'

params {
   genome = 'GRCh37_hmf'
   outdir = '/path/to/output/'
}

process {
   withName: 'REDUX.*' {
      cpus = 32
   }
}
```

The above config items can also be compactly re-written with [dot syntax](https://www.nextflow.io/docs/latest/config.html#blocks) like so:

```
workDir = '/path/to/work/'
params.genome = 'GRCh37_hmf'
params.outdir = '/path/to/output/'
process.withName: 'REDUX.*' { cpus = 32 }
```

### Oncoanalyser arguments as config

The `params` scope is used to define [Oncoanalyser arguments](#oncoanalyser-arguments). Running Oncoanalyser with the 
[above example config](#basic-config-example):
```shell
nextflow run nf-core/oncoanalyser \
-config above_example.config \
# other arguments
```

...is equivalent to running:
```shell
nextflow run nf-core/oncoanalyser \
--genome GRCh37_hmf \
--outdir /path/to/output/ \
# other arguments
```

The `params` scope is also used to define reference data paths (e.g. reference genome, HMFTools resources) as described in 
[Resource files](#resources-files).

### Multiple config files
You may want to keep certain configuration items in separate files. For example:

_oncoanalyser.config_ may contain:
```
params {
   genome = 'GRCh37_hmf'
   outdir = '/path/to/output/'
}
```

...and _processes.config_ may contain:
```
process {
   withName: 'REDUX.*' {
      cpus = 32
   }
}
```

You can provide both when running Oncoanalyser like so:
```shell
nextflow run nf-core/oncoanalyser \
-config oncoanalyser.config \
-config processes.config \
# other arguments
```

## Resource files

### Links

**GRCh37**

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

**GRCh38**

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

### Configuring WGTS resource files

### Configuring panel / targeted sequecing resource files
_TODO_

## Configuring processes

### Compute resources
_TODO_

### Docker images
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

## Acknowledgements

Oncoanalyser was written by Stephen Watts at the [University of Melbourne Centre for Cancer
Research](https://mdhs.unimelb.edu.au/centre-for-cancer-research) with the support of Oliver Hofmann and the Hartwig
Medical Foundation Australia.
