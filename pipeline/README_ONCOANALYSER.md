<h1>
  <img src="nf-core-oncoanalyser_logo_light.png">
</h1>

Oncoanalyser (links: [GitHub](https://github.com/nf-core/oncoanalyser), [nf-core](https://nf-co.re/oncoanalyser/latest/)) is a 
[Nextflow](https://www.nextflow.io/) implementation of the Hartwig Medical Foundation DNA and RNA sequencing analysis pipeline.

**Supported sequencing and sample setups**

| Data type | Sequencing method                                                                | Paired tumor/normal | Tumor-only         |
|:----------|:---------------------------------------------------------------------------------|:--------------------|:-------------------|
| DNA       | Whole genome sequencing (WGS)                                                    | :white_check_mark:  | :white_check_mark: |
| DNA       | Targeted sequencing:<br/> - Whole exome sequencing (WES)<br/> - Panel sequencing | :white_check_mark:  | :white_check_mark: |
| RNA       | Whole transcriptome sequencing (WTS)                                             | -                   | :white_check_mark: |

**Pipeline overview**

The pipeline uses tools from [hmftools](https://github.com/hartwigmedical/hmftools/tree/master/) (except for
[bwa-mem2](https://github.com/bwa-mem2/bwa-mem2), [STAR](https://github.com/alexdobin/STAR) and
[Picard MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard)):

| Task                                                                      | Tool                                                                                                                                                                                                                                                                                                                    |
|:--------------------------------------------------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [Read alignment](#read-alignment)                                         | [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) - DNA<br/>[STAR](https://github.com/alexdobin/STAR) - RNA                                                                                                                                                                                                              |
| [Read post-processing](#read-post-processing)                             | [REDUX](https://github.com/hartwigmedical/hmftools/tree/master/redux) - DNA, duplicate marking and unmapping<br/>[Picard MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard) - RNA, duplicate marking                                                                 |
| [SNV, MNV, INDEL calling](#snv-mnv-indel-calling)                         | [SAGE](https://github.com/hartwigmedical/hmftools/tree/master/sage) - Variant calling<br/>[PAVE](https://github.com/hartwigmedical/hmftools/tree/master/pave) - Transcript/coding effect annotation                                                                                                                     |
| [SV calling](#sv-calling)                                                 | [ESVEE](https://github.com/hartwigmedical/hmftools/tree/master/esvee)<br/>                                                                                                                                                                                                                                              |
| [CNV calling](#cnv-calling)                                               | [AMBER](https://github.com/hartwigmedical/hmftools/tree/master/amber) - B-allele frequencies<br/>[COBALT](https://github.com/hartwigmedical/hmftools/tree/master/cobalt) - Read depth ratios<br/>[PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purple) - Purity/ploidy estimation, variant annotation |
| [SV and driver event interpretation](#sv-and-driver-event-interpretation) | [LINX](https://github.com/hartwigmedical/hmftools/tree/master/linx)                                                                                                                                                                                                                                                     |
| [RNA transcript analysis](#rna-transcript-analysis)                       | [ISOFOX](https://github.com/hartwigmedical/hmftools/tree/master/isofox)                                                                                                                                                                                                                                                 |
| [Oncoviral detection](#oncoviral-detection)                               | [VIRUSbreakend](https://github.com/PapenfussLab/gridss) - Viral content and integration calling<br/>[VirusInterpreter](https://github.com/hartwigmedical/hmftools/tree/master/virus-interpreter) - Post-processing                                                                                                      |
| [Immune analysis](#immune-analysis)                                       | [LILAC](https://github.com/hartwigmedical/hmftools/tree/master/lilac) - HLA typing<br/>[NEO](https://github.com/hartwigmedical/hmftools/tree/master/neo) - Neo-epitope prediction                                                                                                                                       |
| [Mutational signature fitting](#mutational-signature-fitting)             | [SIGS](https://github.com/hartwigmedical/hmftools/tree/master/sigs)                                                                                                                                                                                                                                                     |
| [HRD prediction](#hrd-prediction)                                         | [CHORD](https://github.com/hartwigmedical/hmftools/tree/master/chord)                                                                                                                                                                                                                                                   |
| [Tissue of origin prediction](#tissue-of-origin-prediction)               | [CUPPA](https://github.com/hartwigmedical/hmftools/tree/master/cuppa)                                                                                                                                                                                                                                                   |
| [Summary report](#summary-report)                                         | [ORANGE](https://github.com/hartwigmedical/hmftools/tree/master/orange)                                                                                                                                                                                                                                                 |

## Getting started

This section will assume that: 
- The analysis starts from paired tumor/normal BAMs
- Reads are aligned to the GRCh37 reference genome
- BAMs contain whole genome sequencing data
- Docker images are used to run each tool

> The user has other options including:
> - Starting from FASTQ or other pipeline steps (see: **[Sample sheet](#sample-sheet)**)
> - Using reference genome GRCh38 (see: **[Configuring general resource files](#configuring-general-resource-files)**)
> - Analysing panel sequencing data (see: **[Configuring panel resource files](#configuring-panel-resource-files)**)
> - Using Singularity images (see: **[Container images](#container-images)**)

**1. Install Nextflow**

See: **https://www.nextflow.io/docs/latest/install.html**

**2. Install Docker**

See: **https://docs.docker.com/engine/install/**

**3. Set up resource files**

Download and extract the reference genome and hmftools resources using these **[links](#links)**.

Create a file called `resources.config` which points to the resource file paths:

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
> Jump to: **[Resource files](#resource-files)**, **[Configuration files](#configuration-files)**

**4. Set up sample sheet**

Create a file called `sample_sheet.csv` which points to the sample inputs:

```csv
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
COLO829,COLO829,COLO829T,tumor,dna,bam,/path/to/COLO829T.dna.bam
COLO829,COLO829,COLO829R,normal,dna,bam,/path/to/COLO829R.dna.bam
```

BAM and BAI files for the above COLO829 test sample can be downloaded from [here](test_data/).

> [!TIP]
> Jump to: **[Sample sheet](#sample-sheet)**

**5. Run Oncoanalyser with Nextflow**

```shell
nextflow run nf-core/oncoanalyser \
-profile docker \
-revision pipeline_v6.0 \
-config resources.config \
--mode wgts \
--genome GRCh37_hmf \
--input sample_sheet.csv \
--outdir output/ \
-work-dir output/work
```

> [!TIP]
> Jump to: **[Command line interface](#command-line-interface)**, **[Outputs](#outputs)**, **[Work directory](#work-directory)**

## Table of contents
<!-- TOC -->
  * [Getting started](#getting-started)
  * [Table of contents](#table-of-contents)
  * [Command line interface](#command-line-interface)
    * [Running Oncoanalyser](#running-oncoanalyser)
    * [Nextflow arguments](#nextflow-arguments)
    * [Oncoanalyser arguments](#oncoanalyser-arguments)
  * [Sample sheet](#sample-sheet)
    * [BAM inputs](#bam-inputs)
    * [FASTQ inputs](#fastq-inputs)
    * [Sample modes](#sample-modes)
    * [Multiple patients and/or samples](#multiple-patients-andor-samples)
    * [Running from REDUX BAM](#running-from-redux-bam)
    * [Running specific tools](#running-specific-tools)
  * [Configuration files](#configuration-files)
    * [Basic config example](#basic-config-example)
    * [Oncoanalyser arguments as config](#oncoanalyser-arguments-as-config)
    * [Multiple config files](#multiple-config-files)
  * [Resource files](#resource-files)
    * [Links](#links)
    * [Configuring general resource files](#configuring-general-resource-files)
    * [Configuring panel resource files](#configuring-panel-resource-files)
  * [Configuring processes](#configuring-processes)
    * [Compute resources](#compute-resources)
    * [Maximum resources](#maximum-resources)
    * [Error handling](#error-handling)
  * [Container images](#container-images)
    * [Caching Singularity images](#caching-singularity-images)
    * [Configuring container images](#configuring-container-images)
  * [Outputs](#outputs)
    * [Pipeline information](#pipeline-information)
    * [Read alignment](#read-alignment)
    * [Read post-processing](#read-post-processing)
    * [SNV, MNV, INDEL calling](#snv-mnv-indel-calling)
    * [SV calling](#sv-calling)
    * [CNV calling](#cnv-calling)
    * [SV and driver event interpretation](#sv-and-driver-event-interpretation)
    * [RNA transcript analysis](#rna-transcript-analysis)
    * [Oncoviral detection](#oncoviral-detection)
    * [Immune analysis](#immune-analysis)
    * [Mutational signature fitting](#mutational-signature-fitting)
    * [HRD prediction](#hrd-prediction)
    * [Tissue of origin prediction](#tissue-of-origin-prediction)
    * [Summary report](#summary-report)
  * [Work directory](#work-directory)
  * [Acknowledgements](#acknowledgements)
<!-- TOC -->

## Command line interface

### Running Oncoanalyser
We use the `nextflow run` command to run the Oncoanalyser:

```shell
nextflow run nf-core/oncoanalyser \
-profile docker \
-revision 1.0.0 \
-config hmf_pipeline_resources.config \
--mode wgts \
--genome GRCh37_hmf \
--input sample_sheet.csv \
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

Please section [Outputs](#outputs) and [Work directory](#work-directory) for details on the outputs of Oncoanalyser.

> [!NOTE]
> [Nextflow-specific arguments](https://www.nextflow.io/docs/latest/reference/cli.html) start with a single hyphen (`-`). 
> Oncoanalyser-specific arugments start with two hyphens (`--`).

### Nextflow arguments

All arguments for `nextflow run` are documented in the [CLI reference](https://www.nextflow.io/docs/latest/reference/cli.html#run). The
below table lists some relevant ones:

| Argument&emsp;&emsp; | Description                                                                                                                                                                                                                            |
|:---------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `-config`            | Path to a [configuration file](#configuration-files). [Multiple config files](#multiple-config-files) can be provided.                                                                                                                 |
| `-profile`           | Pre-defined [config profile](https://www.nextflow.io/docs/latest/config.html#config-profiles). For Oncoanalyser, can be `docker`, `singularity`, `test_stub`                                                                           |
| `-latest`            | Pull latest changes before run                                                                                                                                                                                                         |
| `-revision`          | A specific Oncoanalyser branch/tag to run. See the Oncoanalyser [GitHub](https://github.com/nf-core/oncoanalyser) for available branches/tags                                                                                          |
| `-resume`            | [Resume](https://www.nextflow.io/docs/latest/cache-and-resume.html#work-directory) from cached results (by default the previous run). Useful if you've cancelled a run with `CTRL+C`, or a run has crashed and you've fixed the issue. |
| `-stub`              | Dry run. Under the hood, Oncoanalyser runs `touch <outputfile>` rather than actually running the tools. Useful for testing if the arguments and configuration files provided are correct.                                              |
| `-work-dir`          | Path to a directory where Nextflow will put temporary files for each step in the pipeline. If this is not specified, Nextflow will create the `work/` directory in the current directory                                               |
| `-help`              | Show all Nextflow command line arguments and their descriptions                                                                                                                                                                        |

### Oncoanalyser arguments

The below table lists all arguments that can be passed to Oncoanalyser:

| Argument&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; | Description                                                                                                                                                                                                                                                                     |
|:---------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--input`                                                            | Path to a [sample sheet](#sample-sheet)                                                                                                                                                                                                                                         |
| `--outdir`<sup>1</sup>                                               | Path to the output directory. While a process/tool is running, files are temporarily stored in the work directory (see: `-work-dir` [argument](#nextflow-arguments)). Only when the process completes are the files copied to the output directory.                             |
| `--mode`                                                             | Can be:<br/>- `wgts`: Whole genome sequencing and/or whole transcriptome sequencing analysis<br/>- `targeted`: Targeted sequencing analysis (e.g. for panel or whole exome sequencing)                                                                                          |
| `--genome`                                                           | Reference genome version. Can be `GRCh37_hmf` or `GRCh38_hmf`                                                                                                                                                                                                                   |
| `--panel`                                                            | Panel name, e.g. `tso500`                                                                                                                                                                                                                                                       |
| `--force_panel`                                                      | Required flag when `--panel` is not `tso500` (i.e. force run in targeted mode for non-supported panels)                                                                                                                                                                         |
| `--max_cpus`                                                         | Enforce an upper limit of CPUs each process can use, e.g. `16`                                                                                                                                                                                                                  |
| `--max_memory`                                                       | Enforce an upper limit of memory available to each process, e.g. `32.GB`                                                                                                                                                                                                        |
| `--max_time`                                                         | Enforce an upper limit of to the time a process can take, e.g. `240.h`                                                                                                                                                                                                          |
| `--max_fastq_records`                                                | When positive, will use [fastp](https://github.com/OpenGene/fastp) to split fastq files so that each resultant fastq file has no more than max_fastq_records records. When nonpositive, fastp is not used and the provided fastq files are passed as-is to the aligner.         |
| `--processes_exclude`<sup>2</sup>                                    | A comma separated list specifying which processes to skip (e.g. `--processes_exclude lilac,virusinterpreter`). Note: Downstream processes depending on the output of an upstream tool will also be skipped.                                                                     |
| `--processes_include`<sup>2</sup>                                    | When also specifying `--processes_manual`, a comma separated list specifying which processes to include (e.g. `--processes_include lilac,virusinterpreter`). See [Running specific tools](#running-specific-tools) for details on how to set up input files in the sample sheet |
| `--processes_manual`                                                 | Only run processes provided in `--processes_include`                                                                                                                                                                                                                            |
| `--prepare_reference_only`                                           | Only stage reference genome and resource files                                                                                                                                                                                                                                  |
| `--isofox_read_length`                                               | User defined RNA read length used for ISOFOX                                                                                                                                                                                                                                    |
| `--isofox_gc_ratios`                                                 | User defined ISOFOX expected GC ratios file                                                                                                                                                                                                                                     |
| `--isofox_counts`                                                    | User defined ISOFOX expected counts files (read length dependent)                                                                                                                                                                                                               |
| `--isofox_tpm_norm`                                                  | User defined ISOFOX TPM normalisation file for panel data                                                                                                                                                                                                                       |
| `--isofox_gene_ids`                                                  | User defined ISOFOX gene list file for panel data.                                                                                                                                                                                                                              |
| `--isofox_functions`                                                 | Semicolon-separated list of ISOFOX functions to run. Default: `TRANSCRIPT_COUNTS;ALT_SPLICE_JUNCTIONS;FUSIONS;RETAINED_INTRONS`                                                                                                                                                 |
| `--fastp_umi`                                                        | Enable UMI processing by fastp                                                                                                                                                                                                                                                  |
| `--fastp_umi_location`                                               | Passed to fastp arg `--umi_loc`. Can be `per_index` or `per_read`                                                                                                                                                                                                               |
| `--fastp_umi_length`                                                 | Passed to fastp arg `--umi_len`. Expected length (number of bases) of the UMI                                                                                                                                                                                                   |
| `--fastp_umi_skip`                                                   | Passed to fastp arg `--umi_skip`. Number of bases to skip following UMI                                                                                                                                                                                                         |
| `--redux_umi`                                                        | Enable UMI processing by REDUX                                                                                                                                                                                                                                                  |
| `--redux_umi_duplex_delim`                                           | UMI duplex delimiter as used by REDUX, Default: `_`                                                                                                                                                                                                                             |
| `--ref_data_hmf_data_path`                                           | Path to hmftools resource files                                                                                                                                                                                                                                                 |
| `--ref_data_panel_data_path`                                         | Path to panel resource files                                                                                                                                                                                                                                                    |
| `--ref_data_hla_slice_bed`                                           | Path to HLA slice BED file                                                                                                                                                                                                                                                      |
| `--create_stub_placeholders`                                         | Create placeholders for resource files during stub run                                                                                                                                                                                                                          |
| `--email`                                                            | Email address for completion summary                                                                                                                                                                                                                                            |
| `--monochrome_logs`                                                  | Do not use coloured log outputs                                                                                                                                                                                                                                                 |

Notes:
1. **WARNING**: Cannot be provided to a [config file](#configuration-files)
2. Valid process names are: `alignment`, `amber`, `bamtools`, `chord`, `cobalt`, `cuppa`, `esvee`, `isofox`, `lilac`, `linx`, `neo`,
   `orange`, `pave`, `purple`, `redux`, `sage`, `sigs`, `virusinterpreter`

## Sample sheet

The sample sheet is a comma separated table with the following columns:
- `subject_id`: Top level grouping
- `group_id`: Groups `sample_id` entries (e.g. group tumor DNA, normal DNA, tumor RNA) into the same analysis 
- `sample_id`
- `sample_type`: `tumor` or `normal`
- `sequence_type`: `dna` or `rna`
- `filetype`: `bam`, `bai`, `fastq`, or see **[Running specific tools](#running-specific-tools)** for other valid values
- `filepath`: Absolute filepath to input file. Can be local filepath, URL, or S3 URI
- `info`: Sequencing library and lane info for **[FASTQ inputs](#fastq-inputs)**

### BAM inputs
Below is an example sample sheet with BAM files for a tumor/normal WGS run:

```csv
subject_id,group_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.bam
PATIENT1,PATIENT1,PATIENT1-R,normal,dna,bam,/path/to/PATIENT1-R.dna.bam
```

BAM indexes (.bai files) are expected to be in the same directory as the BAM files. Alternatively, provide the BAM index path by 
providing `bai` under column `filetype`:

```csv
subject_id,group_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.bam
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bai,/path/to/PATIENT1-T.dna.bam.bai
PATIENT1,PATIENT1,PATIENT1-R,normal,dna,bam,/path/to/PATIENT1-R.dna.bam
PATIENT1,PATIENT1,PATIENT1-R,normal,dna,bai,/path/to/PATIENT1-R.dna.bam.bai
```

### FASTQ inputs
Below is an example sample sheet with FASTQ files for a tumor/normal WGS run:

```csv
subject_id,group_id,sample_id,sample_type,sequence_type,filetype,filepath,info
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,fastq,/path/to/PATIENT1-T_S1_L001_R1_001.fastq.gz;/path/to/PATIENT1-T_S1_L001_R2_001.fastq.gz,library_id:S1;lane:001
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,fastq,/path/to/PATIENT1-T_S1_L002_R1_001.fastq.gz;/path/to/PATIENT1-T_S1_L002_R2_001.fastq.gz,library_id:S1;lane:002
PATIENT1,PATIENT1,PATIENT1-R,normal,dna,fastq,/path/to/PATIENT1-R_S2_L001_R1_001.fastq.gz;/path/to/PATIENT1-R_S2_L001_R2_001.fastq.gz,library_id:S2;lane:001
PATIENT1,PATIENT1,PATIENT1-R,normal,dna,fastq,/path/to/PATIENT1-R_S2_L002_R1_002.fastq.gz;/path/to/PATIENT1-R_S2_L002_R2_001.fastq.gz,library_id:S2;lane:002
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
subject_id,group_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.bam
```

**Tumor-only RNA**
```csv
subject_id,group_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T-RNA,tumor,rna,bam,/path/to/PATIENT1-T.rna.bam
```

**Tumor/normal DNA, tumor-only RNA**

```csv
subject_id,group_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.bam
PATIENT1,PATIENT1,PATIENT1-R,normal,dna,bam,/path/to/PATIENT1-R.dna.bam
PATIENT1,PATIENT1,PATIENT1-T-RNA,tumor,dna,bam,/path/to/PATIENT1-T.rna.bam
```

### Multiple patients and/or samples

Suppose you have multiple patients, each with one or more biopsies taken from different years.

You could then set:
- `subject_id` to the patient ID
- `group_id` to the set of samples for a particular year (e.g. `PATIENT1-YEAR1`)
- `sample_id` to the actual sample IDs in the sample set for that year

For example:

```csv
subject_id,group_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1-YEAR1,PATIENT1-YEAR1-T,tumor,dna,bam,/path/to/PATIENT1-YEAR1-T.dna.bam
PATIENT1,PATIENT1-YEAR1,PATIENT1-YEAR1-R,normal,dna,bam,/path/to/PATIENT1-YEAR1-R.dna.bam
PATIENT1,PATIENT1-YEAR2,PATIENT1-YEAR2-T,tumor,dna,bam,/path/to/PATIENT1-YEAR2-T.dna.bam
PATIENT1,PATIENT1-YEAR2,PATIENT1-YEAR2-R,normal,dna,bam,/path/to/PATIENT1-YEAR2-R.dna.bam
PATIENT2,PATIENT2-YEAR1,PATIENT2-YEAR1-T,tumor,dna,bam,/path/to/PATIENT2-YEAR1-T.dna.bam
PATIENT2,PATIENT2-YEAR1,PATIENT2-YEAR1-R,normal,dna,bam,/path/to/PATIENT2-YEAR1-R.dna.bam
```

### Running from REDUX BAM
For DNA sequencing analyses, read alignment with [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) and read pre-processing with 
[REDUX](https://github.com/hartwigmedical/hmftools/tree/master/redux) are the pipeline steps that take the most time and compute resources. 
Thus, we can re-run Oncoanalyser from a REDUX BAM if it is already exists, e.g. due to updates to downstream tools from 
[hmftools](https://github.com/hartwigmedical/hmftools/tree/master/).

Simply provide the REDUX BAM path, specifying `bam_redux` under `filetype`:

```csv
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam_redux,/path/to/PATIENT1-T.dna.redux.bam
```

The `*.jitter_params.tsv` and `*.ms_table.tsv.gz` REDUX output files are expected to be in the same directory as the REDUX BAM. If these 
files are located elsewhere, their paths can also be explicitly provided by specifying `redux_jitter_tsv` and `redux_ms_tsv` under `filetype`:

```csv
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam_redux,/path/to/PATIENT1-T.dna.redux.bam
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,redux_jitter_tsv,/path/to/PATIENT1-T.dna.jitter_params.tsv
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,redux_ms_tsv,/path/to/PATIENT1-T.dna.ms_table.tsv.gz
```

### Running specific tools
It is possible to run Oncoanalyser from any tool from [hmftools](https://github.com/hartwigmedical/hmftools/tree/master/). For example, you may want to 
run [CUPPA](https://github.com/hartwigmedical/hmftools/tree/master/cuppa) and already have the outputs from
[PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purple), 
[LINX](https://github.com/hartwigmedical/hmftools/tree/master/linx), 
and [VirusInterpreter](https://github.com/hartwigmedical/hmftools/tree/master/virus-interpreter). In this case, you would provide the 
outputs from those tools to the sample sheet, specifying entries where `filetype` is `purple_dir`, `linx_anno_dir`, and 
`virusinterpreter_dir`:

```csv
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,purple_dir,/path/to/purple/dir/
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,linx_anno_dir,/path/to/linx/dir/
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,virusinterpreter_dir,/path/to/virus/dir/
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
- [General resource files](#configuring-general-resource-files)
- [Panel resource files](#configuring-panel-resource-files)
- [Compute resources](#compute-resources)
- [Maximum resources](#maximum-resources)
- [Error handling](#error-handling)
- [Container images](#configuring-container-images)

> [!NOTE]
> Configuration is fully detailed in the **[Nextflow](https://www.nextflow.io/docs/latest/config.html)** and 
> **[nf-core](https://nf-co.re/docs/usage/getting_started/configuration)** documentation

### Basic config example
Config items can be declared using [blocks](#https://www.nextflow.io/docs/latest/config.html#blocks), where curly brackets define the scope 
of the encapsulated config items. The below example has the `params` scope, with `workDir` being 
[un-scoped](https://www.nextflow.io/docs/latest/reference/config.html#unscoped-options):

```
params {
   ref_data_hmf_data_path = '/path/to/hmf_pipeline_resources/'
   redux_umi = true
}

workDir = '/path/to/work/'
```

The above config items can also be compactly re-written with [dot syntax](https://www.nextflow.io/docs/latest/config.html#blocks) like so:

```
params.ref_data_hmf_data_path = '/path/to/hmf_pipeline_resources/'
params.redux_umi = true

workDir = '/path/to/work/'
```

### Oncoanalyser arguments as config

The `params` scope can be used to define [Oncoanalyser arguments](#oncoanalyser-arguments). Running Oncoanalyser with the 
[above example config](#basic-config-example):
```shell
nextflow run nf-core/oncoanalyser \
-config above_example.config \
# other arguments
```

...is equivalent to running:
```shell
nextflow run nf-core/oncoanalyser \
--ref_data_hmf_data_path /path/to/hmf_pipeline_resources/ \
--redux_umi \
# other arguments
```

The `params` scope is also used to define reference data paths (e.g. reference genome, hmftools resources) as detailed in 
[Resource files](#resource-files).

### Multiple config files
You may want to keep certain configuration items in separate files. For example:

_resource_files.config_ may contain:
```
params {
   ref_data_hmf_data_path = '/path/to/hmf_pipeline_resources/'
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
-config resource_files.config \
-config processes.config \
# other arguments
```

## Resource files

### Links

**GRCh37**

| Type         | Description          | Name                                                                                                                                                                                                |
|--------------|----------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| hmftools     | hmftools resources   | [hmf_pipeline_resources.37_v6.0--2.tar.gz](https://storage.googleapis.com/hmf-public/HMFtools-Resources/oncoanalyser/v6_0/37/hmf_pipeline_resources.37_v6.0--2.tar.gz)                              |
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
| hmftools     | hmftools resources   | [hmf_pipeline_resources.38_v6.0--2.tar.gz](https://storage.googleapis.com/hmf-public/HMFtools-Resources/oncoanalyser/v6_0/38/hmf_pipeline_resources.38_v6.0--2.tar.gz)                                              |
| Genome       | FASTA                | [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.0/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna)                                      |
| Genome       | FASTA index          | [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.0/samtools_index/1.16/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai)          |
| Genome       | FASTA seq dictionary | [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.dict](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.0/samtools_index/1.16/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.dict)        |
| Genome       | bwa-mem2 index image | [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.img](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.0/bwa_index_image/0.7.17-r1188/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.img) |
| Genome       | bwa-mem2 index       | [bwa-mem2_index/2.2.1.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.1/bwa-mem2_index/2.2.1.tar.gz)                                                                              |
| Genome       | GRIDSS index         | [gridss_index/2.13.2.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.1/gridss_index/2.13.2.tar.gz)                                                                                |
| Genome (RNA) | STAR index           | [star_index/gencode_38/2.7.3a.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.0/star_index/gencode_38/2.7.3a.tar.gz)                                                              |
| Panel        | TSO500 data          | [panels/tso500_5.34_38--1.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/panels/tso500_5.34_38--1.tar.gz)                                                                           |

### Configuring general resource files

The below example shows the most essential config items when configuring resource files. Not all items are required depending on the 
experimental setup. Please see the inline comments for details.

> [!NOTE]
> Single line comments start with `//`. Multi-line comments start with `/*` and end with `*/`

```
params {
   genomes {
   
      'GRCh37_hmf' { // Can be 'GRCh37_hmf' or 'GRCh38_hmf'
      
         fasta         = "/path/to/Homo_sapiens.GRCh37.GATK.illumina.fasta"
         fai           = "/path/to/Homo_sapiens.GRCh37.GATK.illumina.fasta.fai"
         dict          = "/path/to/Homo_sapiens.GRCh37.GATK.illumina.fasta.dict"
         img           = "/path/to/Homo_sapiens.GRCh37.GATK.illumina.fasta.img"
         
         // Required if aligning reads from FASTQ files (can be skipped when running from BAM files)
         bwamem2_index = "/path/to/bwa-mem2_index/"
         
         // Required if running VIRUSbreakend
         gridss_index  = "/path/to/gridss_index/"
         
         // Required only for RNA sequencing data
         star_index    = "/path/to/star_index/"
      }
      
      // Both GRCh37_hmf and GRCh38_hmf entries can be provided!
      'GRCh38_hmf' {
         fasta         = "/path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
         // Provide remaining options in a similar manner as for 'GRCh37_hmf' above
      }
      
   }
   
   // Always required
   ref_data_hmf_data_path = "/path/to/hmf_pipeline_resources/"
}
```

### Configuring panel resource files
Running in `--mode targeted`, requires some additional resources files to be configured:

```
params {   
   ref_data_panel_data_path = "/path/to/panel_resources/"
   
   // These are relative paths within the dir provided by `ref_data_panel_data_path`
   panel_data_paths {
      
      custom_panel { // This is the name that should be passed to the `--panel` argument
         
         // Can be '37' or '38'
         '37' {
              
            driver_gene_panel           = 'common/custom_panel.driver_gene_panel.tsv'
            sage_actionable_panel       = 'variants/custom_panel.coding_panel.v37.bed.gz'
            sage_coverage_panel         = 'variants/custom_panel.coverage_panel.v37.bed.gz'
            pon_artefacts               = 'variants/custom_panel.sage_pon_artefacts.tsv.gz'
            target_region_bed           = 'custom_panel.panel_regions.v37.bed.gz'
            target_region_normalisation = 'copy_number/custom_panel.cobalt_normalisation.37.tsv'
            target_region_ratios        = 'copy_number/custom_panel.target_regions_ratios.37.tsv'
            target_region_msi_indels    = 'copy_number/custom_panel.target_regions_msi_indels.37.tsv'
            
            // Optional. These can be omitted, or provided a falsy value such as '' or []
            isofox_tpm_norm             = ''
            isofox_gene_ids             = ''
            isofox_counts               = ''
            isofox_gc_ratios            = ''
         }
      }
   }
}
```

When running Oncoanalyser: 
- Provide both the general and panel resources config files to `-config`
- Pass the panel name to `--panel`. This should match the name defined in the panel resources config file
- Provide argument `--force_panel` if `--panel` is not `tso500` (this is currently the only supported panel type) 

```shell
nextflow run nf-core/oncoanalyser \
--panel custom_panel \
--force_panel \
-config general_resources.config \
-config panel_resources.config \
--mode targeted \
# other arguments
```

## Configuring processes

There are many options for configuration processes. However, this section will go over some common use cases.

> [!NOTE]
> Configuration of processes is fully detailed in the [nf-core](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources).
> All options (a.k.a 'directives') for configuring processes are detailed in the [Nextflow process reference docs](https://www.nextflow.io/docs/latest/reference/process.html#directives).

### Compute resources

Each hmftool is run within a Nextflow **process**. We can use the `process` scope and `withName` to select tools by name and set compute 
resources options (as well as other [config options](https://www.nextflow.io/docs/latest/reference/process.html)):

```
process {
   withName: 'SAGE_SOMATIC' {
      cpus = 32
      memory = 128.GB
      disk = 1024.GB
      time = 48.h
   }
}
```

Values with a units are provided in quotes with a space or without quotes using a dot, e.g. `'128 GB'` or `128.GB`.

Please see the [Nextflow process reference docs](https://www.nextflow.io/docs/latest/reference/process.html#directives) to see all possible
options. The following links is the documentation for the ones used above:
**[cpus](https://www.nextflow.io/docs/latest/reference/process.html#cpus)**,
**[memory](https://www.nextflow.io/docs/latest/reference/process.html#memory)**,
**[time](https://www.nextflow.io/docs/latest/reference/process.html#time)**,
**[disk](https://www.nextflow.io/docs/latest/reference/process.html#disk)**.

We can also use a regular expression to select multiple processes. SAGE for example has the processes `SAGE_SOMATIC`, `SAGE_GERMLINE` and
`SAGE_APPEND`. We can select all 3 like so:

```
process {
   withName: 'SAGE.*' {
      cpus = 32
   }
}
```

Processes are also grouped by compute resource **labels**, with the main ones being (in order of increasing compute load) `process_single`,
`process_low`, `process_medium` and `process_high`. The labels `process_medium_memory` and `process_high_memory` are only used for creating
genome indexes. We can use `withLabel` to set options for all tools with an associated label:

```
process {
   withLabel: 'process_low' {
      cpus = 2
   }
}
```

### Maximum resources

The [maximum resources](https://www.nextflow.io/docs/latest/reference/process.html#resourcelimits) for any process can also be set using 
`resourceLimits`. If a process requests more resources than allowed (e.g. a process requests 64 cores but the largest node in a cluster has 32), 
the process would normally fail or cause the pipeline to hang forever as it will never be scheduled. Setting `resourceLimits` will 
automatically reduce the process resources to comply with the provided limits before the job is submitted.
```
process {
   resourceLimits = [
      cpus: 32,
      memory: 128.GB,
      time: 48.h
   ]
}
```

### Error handling

We can use `errorStrategy` and `maxRetries` to determine how Oncoanalyser proceeds when encountering an error. For example, to retry 3 times 
on any error for any process, we can provide this config:

```
process {
   errorStrategy = 'retry'
   maxRetries = 3
}
```

Valid values for `errorStrategy` are (details in the Nextflow [documentation](https://www.nextflow.io/docs/latest/reference/process.html#errorstrategy)):
- `retry`: Retry the process
- `terminate`: Fail the pipeline immediately
- `finish`: Terminate after submitted and running processes are done
- `ignore`: Allow the pipeline to continue upon error

Process selectors can also be used to target specific processes for error handling:
```
process {
   withName: 'SAGE_SOMATIC' {
      errorStrategy = 'retry'
      maxRetries = 3
   }
}
```

## Container images

Oncoanalyser by default uses **[Docker](https://www.docker.com/)** and [**Singularity**](https://docs.sylabs.io/guides/3.0/user-guide/quick_start.html#) 
images built by the **[bioconda-recipes](https://github.com/bioconda/bioconda-recipes/tree/master/recipes)** Azure CI/CD infrastructure.

Use `-profile docker` or `-profile singularity` to tell Oncoanalyser whether to run with Docker or Singularity respectively. For example:
```shell
nextflow run nf-core/oncoanalyser 
-profile docker \
# other arguments
```

Docker images built by Hartwig's Google Cloud CI/CD infrastructure are also available (though not used by default by Oncoanalyser).

> [!NOTE]
> Configuration of container images is fully detailed in the **[nf-core](https://nf-co.re/docs/usage/configuration#docker-registries)** and 
> **[Nextflow](https://www.nextflow.io/docs/edge/container.html#containers)** documentation

Docker and singularity image URIs/URLs have consistent patterns:

| Source   | Platform    | Host                                                           | URI or URL                                                                                                                                                                         |
|:---------|:------------|:---------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Bioconda | Docker      | [quay.io](https://quay.io/organization/biocontainers)          | **Pattern**: `quay.io/biocontainers/hmftools-{TOOL}:{TAG}`<br/> **Example**: `quay.io/biocontainers/hmftools-sage:4.0_beta--hdfd78af_4`                                            |
| Bioconda | Singularity | [Galaxy Project](https://depot.galaxyproject.org/singularity/) | **Pattern**: `https://depot.galaxyproject.org/singularity/hmftools-{tool}:{tag}` <br/> **Example**: https://depot.galaxyproject.org/singularity/hmftools-sage:4.0_beta--hdfd78af_4 |
| Hartwig  | Docker      | [Dockerhub](https://hub.docker.com/r/hartwigmedicalfoundation) | **Pattern**: `docker.io/hartwigmedicalfoundation/{TOOL}:{TAG}` <br/> **Example**: `docker.io/hartwigmedicalfoundation/sage:4.0-rc.2`                                               |

Bioconda recipes also have a consistent URL pattern:
- **Pattern**: `https://github.com/bioconda/bioconda-recipes/tree/master/recipes/hmftools-{tool}`
- **Example**: https://github.com/bioconda/bioconda-recipes/tree/master/recipes/hmftools-sage

These patterns are useful to know as the [bioconda-recipes](https://github.com/bioconda/bioconda-recipes/tree/master/recipes),
[quay.io](https://quay.io/organization/biocontainers), and [Galaxy Project](https://depot.galaxyproject.org/singularity/) repos especially 
have thousands of entries but poor search functionality.

Oncoanalyser/Nextflow automatically pulls Bioconda images at runtime. However, images can also be manually pulled from URIs/URLs. For example:
```shell
## Docker: Downloads into your local Docker repository 
docker pull quay.io/biocontainers/hmftools-sage:4.0_beta--hdfd78af_4

## Singularity: Downloads image to a file called 'hmftools-sage:4.0_beta--hdfd78af_4' 
singularity pull https://depot.galaxyproject.org/singularity/hmftools-sage:4.0_beta--hdfd78af_4
```

### Caching Singularity images

Some compute environments, especially HPCs (high performance clusters), grant limited network access which prevents Oncoanalyser/Nextflow 
from automatically pulling images at runtime. To get around this, we can manually download the Singularity images to a directory: 

```shell
cd /path/to/cache_dir/

## For the image name provided to the `--name` argument, remove 'https://' and replace '/' with '-'
singularity pull \
--name depot.galaxyproject.org-singularity-hmftools-sage:4.0_beta--hdfd78af_4.img \
https://depot.galaxyproject.org/singularity/hmftools-sage:4.0_beta--hdfd78af_4

## Repeat for all singularity images that Oncoanalyser uses
## singularity pull ...
```

and set the `NXF_SINGULARITY_CACHEDIR` environment variable ([Nextflow documentation](https://www.nextflow.io/docs/edge/container.html#singularity-docker-hub))
to tell Oncoanalyser/Nextflow where to look for local images at runtime:

```shell
export NXF_SINGULARITY_CACHEDIR='/path/to/cache_dir/'

nextflow run nf-core/oncoanalyser 
-profile singularity \
# ...other arguments
```

Alternatively, the path to the Singularity cache dir can also be provided to a config file:

```
singularity {
    cacheDir = '/path/to/cache_dir/'
}
```

to be passed to `nextflow run`:
```shell
export NXF_SINGULARITY_CACHEDIR='/path/to/cache_dir/'

nextflow run nf-core/oncoanalyser 
-profile singularity \
-config singularity.config \
# ...other arguments
```

> [!NOTE]
> All singularity options are detailed in the **[Singularity Nextflow documentation](https://www.nextflow.io/docs/latest/reference/config.html#config-singularity)**

### Configuring container images

We can override the default container image used by Oncoanalyser like so:

```
process {
   withName: 'SAGE.*' {
      container = docker.io/hartwigmedicalfoundation/sage:4.0-rc.2
   }
   
  withName: 'ESVEE.*' {
     container = docker.io/hartwigmedicalfoundation/esvee:1.0-rc.4
  }
}
```

This is useful for example when you want to use updated container images that are not yet officially supported (e.g. betas or release candidates). 

In general, the process names for all hmftools are `{TOOL}` or `{TOOL}_{SUBPROCESS}`. For example, SAGE has the processes: `SAGE_SOMATIC`, 
`SAGE_GERMLINE`, `SAGE_APPEND`. Therefore, use regex suffix `.*` (e.g. `SAGE.*`) to capture the subprocesses for each tool.

## Outputs

Oncoanalyser writes output files to the below directory tree structure at the path provided by the `--outdir` 
[argument](#oncoanalyser-arguments). Files are grouped by the `group_id` provided in the [sample sheet](#sample-sheet), then by tool:

```shell
output/
├── pipeline_info/
├── group_id_1/
│   ├── alignments/
│   ├── amber/
│   ├── bamtools/
│   ├── chord/
│   ├── cobalt/
│   ├── cuppa/
│   ├── esvee/
│   ├── isofox/
│   ├── lilac/
│   ├── linx/
│   ├── orange/
│   ├── pave/
│   ├── purple/
│   ├── sage/
│   ├── sigs/
│   ├── virusbreakend/
│   └── virusinterpreter/
│   
├── group_id_2/
│   └── ...
│   
...
```

All intermediate files used by each process are kept in the Nextflow [work directory](#work-directory). Once an analysis has completed this 
directory can be removed.

### Pipeline information

Created by Nextflow:

```shell
pipeline_info/
├── execution_report_<date_time>.html   # HTML report of execution metrics and details
├── execution_timeline_<date_time>.html # Timeline diagram showing process start/duration/finish
├── execution_trace_<date_time>.txt     # Resource usage
├── pipeline_dag_<date_time>.html       # Pipeline diagram showing how each process is connected
```

Created by Oncoanalyser:
```shell
├── params_<date_time>.json             # Parameters used by the pipeline run
└── software_versions.yml               # Tool versions
```

### Read alignment

No outputs from **bwa-mem2** and **STAR** are published.

### Read post-processing

**REDUX**: Duplicate marking and unmapping

```shell
<group_id>/alignments/
├── dna
│   ├── <tumor_dna_id>.jitter_params.tsv         # Microsatellite jitter model parameters
│   ├── <tumor_dna_id>.ms_table.tsv.gz           # Aggregated repeat units and repeat counts
│   ├── <tumor_dna_id>.redux.bam                 # Read alignments
│   ├── <tumor_dna_id>.redux.bam.bai             # Read alignments index
│   ├── <tumor_dna_id>.redux.duplicate_freq.tsv  # Duplicate read frequencies
│   ├── <tumor_dna_id>.repeat.tsv.gz             # Repeat units and repeat counts per site
│   ├── <normal_dna_id>.jitter_params.tsv        # See above
│   ├── <normal_dna_id>.ms_table.tsv.gz          # See above
│   ├── <normal_dna_id>.redux.bam                # See above
│   ├── <normal_dna_id>.redux.bam.bai            # See above
│   ├── <normal_dna_id>.redux.duplicate_freq.tsv # See above
│   └── <normal_dna_id>.repeat.tsv.gz            # See above
```

**Picard MarkDuplicates**: Duplicate marking

```shell
└── rna
    ├── `<tumor_rna_id>.md.bam`     # Read alignments               
    ├── `<tumor_rna_id>.md.bam.bai` # Read alignments index         
    └── `<tumor_rna_id>.md.metrics` # Duplicate read marking metrics
```

### SNV, MNV, INDEL calling

**SAGE**: Variant calling

```shell
<group_id>/sage
├── somatic
│   ├── <normal_dna_id>.sage.bqr.png            # Normal DNA sample base quality recalibration metrics plot
│   ├── <normal_dna_id>.sage.bqr.tsv            # Normal DNA sample base quality recalibration metrics
│   ├── <tumor_dna_id>.sage.bqr.png             # Tumor DNA sample base quality recalibration metrics plot
│   ├── <tumor_dna_id>.sage.bqr.tsv             # Tumor DNA sample base quality recalibration metrics
│   ├── <tumor_dna_id>.sage.exon.medians.tsv    # Tumor DNA sample exon median depths
│   ├── <tumor_dna_id>.sage.gene.coverage.tsv   # Tumor DNA sample gene coverages
│   ├── <tumor_dna_id>.sage.somatic.vcf.gz      # Tumor DNA sample small variant calls
│   └── <tumor_dna_id>.sage.somatic.vcf.gz.tbi  # Tumor DNA sample small variant calls index
├── germline
│   ├── <normal_dna_id>.sage.bqr.png            # Tumor DNA sample base quality recalibration metrics plot             
│   ├── <normal_dna_id>.sage.bqr.tsv            # Tumor DNA sample base quality recalibration metrics
│   ├── <normal_dna_id>.sage.exon.medians.tsv   # Normal DNA sample exon median depths
│   ├── <normal_dna_id>.sage.gene.coverage.tsv  # Normal DNA sample gene coverages
│   ├── <tumor_dna_id>.sage.bqr.png             # Normal DNA sample base quality recalibration metrics plot
│   ├── <tumor_dna_id>.sage.bqr.tsv             # Normal DNA sample base quality recalibration metrics
│   ├── <tumor_dna_id>.sage.germline.vcf.gz     # Normal DNA sample filtered small variant calls
│   └── <tumor_dna_id>.sage.germline.vcf.gz.tbi # Normal DNA sample filtered small variant calls index
└── append
    ├── <normal_dna_id>.sage.append.vcf.gz      # Normal VCF with SMNVs and RNA data appended
    └── <tumor_dna_id>.sage.append.vcf.gz       # Tumor VCF with SMNVs and RNA data appended
```

**PAVE**: Transcript/coding effect annotation

```shell
<group_id>/pave/
├── <tumor_dna_id>.sage.germline.pave.vcf.gz     # VCF with annotated germline SAGE SMNVs
├── <tumor_dna_id>.sage.germline.pave.vcf.gz.tbi # VCF index
├── <tumor_dna_id>.sage.somatic.pave.vcf.gz      # VCF with annotated somatic SAGE SMNVs
└── <tumor_dna_id>.sage.somatic.pave.vcf.gz.tbi  # VCF index
```

### SV calling

**ESVEE**: Variant calling

```shell
<group_id>/esvee/
├── prep
│   ├── <tumor_dna_id>.esvee.prep.bam                 # BAM with candidate SV reads 
│   ├── <tumor_dna_id>.esvee.prep.bam.bai             # BAM index
│   ├── <tumor_dna_id>.esvee.prep.disc_stats.tsv      # Discordant reads stats
│   ├── <tumor_dna_id>.esvee.prep.fragment_length.tsv # Fragment length stats
│   ├── <tumor_dna_id>.esvee.prep.junction.tsv        # Candidate junctions
│   ├── <normal_dna_id>.esvee.prep.bam                # BAM with candidate SV reads
│   └── <normal_dna_id>.esvee.prep.bam.bai            # BAM index
├── assemble
│   ├── <tumor_dna_id>.esvee.assembly.tsv             # Breakend assemblies
│   ├── <tumor_dna_id>.esvee.alignment.tsv            # Assemblies realigned to the ref genome
│   ├── <tumor_dna_id>.esvee.breakend.tsv             #
│   ├── <tumor_dna_id>.esvee.phased_assembly.tsv      #
│   ├── <tumor_dna_id>.esvee.raw.vcf.gz               # VCF with candidate breakends
│   └── <tumor_dna_id>.esvee.raw.vcf.gz.tbi           # VCF with candidate breakends
├── depth_annotation
│   ├── <tumor_dna_id>.esvee.ref_depth.vcf.gz         # VCF annotated with normal sample read depths
│   └── <tumor_dna_id>.esvee.ref_depth.vcf.gz.tbi     # VCF index
└── caller
    ├── <tumor_dna_id>.esvee.germline.vcf.gz          # VCF with germline breakends
    ├── <tumor_dna_id>.esvee.germline.vcf.gz.tbi      # VCF index
    ├── <tumor_dna_id>.esvee.somatic.vcf.gz           # VCF with somatic breakends
    ├── <tumor_dna_id>.esvee.somatic.vcf.gz.tbi       # VCF index
    ├── <tumor_dna_id>.esvee.unfiltered.vcf.gz        # VCF with unfiltered breakends
    └── <tumor_dna_id>.esvee.unfiltered.vcf.gz.tbi    # VCF index
```

### CNV calling

**AMBER**: B-allele frequencies

```shell
<group_id>/amber/
├── <tumor_dna_id>.amber.baf.pcf                  # Piecewise constant fit on B-allele frequencies
├── <tumor_dna_id>.amber.baf.tsv.gz               # B-allele frequencies
├── <tumor_dna_id>.amber.contamination.tsv        # Contamination TSV
├── <tumor_dna_id>.amber.contamination.vcf.gz     # Contamination sites
├── <tumor_dna_id>.amber.contamination.vcf.gz.tbi # Sample contamination sites index
├── <tumor_dna_id>.amber.qc                       # QC file
├── <normal_dna_id>.amber.homozygousregion.tsv    # Regions of homozygosity
├── <normal_dna_id>.amber.snp.vcf.gz              # SNP sites VCF
├── <normal_dna_id>.amber.snp.vcf.gz.tbi          # VCF index
└── amber.version                                 # Tool version
```

**COBALT**: Read depth ratios

```shell
<group_id>/cobalt/
├── <tumor_dna_id>.cobalt.gc.median.tsv     # GC median read depths
├── <tumor_dna_id>.cobalt.ratio.pcf         # Piecewise constant fit
├── <tumor_dna_id>.cobalt.ratio.tsv.gz      # Read counts and ratios (with reference or supposed diploid)
├── <normal_dna_id>.cobalt.gc.median.tsv    # GC median read depths
├── <normal_dna_id>.cobalt.ratio.median.tsv # Chromosome median ratios  
├── <normal_dna_id>.cobalt.ratio.pcf        # Piecewise constant fit
└── cobalt.version                          # Tool version
```

**PURPLE**: Purity/ploidy estimation, variant annotation

```shell
<group_id>/purple/
├── <tumor_dna_id>.purple.cnv.gene.tsv                # Somatic gene copy number
├── <tumor_dna_id>.purple.cnv.somatic.tsv             # Copy number variant segments
├── <tumor_dna_id>.purple.driver.catalog.germline.tsv # Germline DNA sample driver events
├── <tumor_dna_id>.purple.driver.catalog.somatic.tsv  # Somatic DNA sample driver events
├── <tumor_dna_id>.purple.germline.deletion.tsv       # Germline DNA deletions
├── <tumor_dna_id>.purple.germline.vcf.gz             # Germline SAGE SMNVs with PURPLE annotations
├── <tumor_dna_id>.purple.germline.vcf.gz.tbi         # VCF index
├── <tumor_dna_id>.purple.purity.range.tsv            # Purity/ploidy model fit scores across a range of purity values
├── <tumor_dna_id>.purple.purity.tsv                  # Purity/ploidy summary
├── <tumor_dna_id>.purple.qc                          # QC file
├── <tumor_dna_id>.purple.segment.tsv                 # Genomic copy number segments
├── <tumor_dna_id>.purple.somatic.clonality.tsv       # Clonality peak model data
├── <tumor_dna_id>.purple.somatic.hist.tsv            # Somatic variants histogram data
├── <tumor_dna_id>.purple.somatic.vcf.gz              # Tumor SAGE SMNVs with PURPLE annotations
├── <tumor_dna_id>.purple.somatic.vcf.gz.tbi          # VCF index
├── <tumor_dna_id>.purple.sv.germline.vcf.gz          # Germline ESVEE SVs with PURPLE annotations
├── <tumor_dna_id>.purple.sv.germline.vcf.gz.tbi      # VCF index
├── <tumor_dna_id>.purple.sv.vcf.gz                   # Somatic ESVEE SVs with PURPLE annotations
├── <tumor_dna_id>.purple.sv.vcf.gz.tbi               # VCF index
├── circos/         # Circos plot data
├── plot/           # PURPLE plots
└── purple.version  # Tool version
```

### SV and driver event interpretation

**LINX**: SV and driver event interpretation

```shell
<group_id>/linx/
├── germline_annotations
│   ├── <tumor_dna_id>.linx.germline.breakend.tsv       # Normal sample breakend data
│   ├── <tumor_dna_id>.linx.germline.clusters.tsv       # Normal sample clustered events
│   ├── <tumor_dna_id>.linx.germline.disruption.tsv     # 
│   ├── <tumor_dna_id>.linx.germline.driver.catalog.tsv # Normal sample driver events
│   ├── <tumor_dna_id>.linx.germline.links.tsv          # 
│   ├── <tumor_dna_id>.linx.germline.svs.tsv            #
│   └── linx.version                                    # Tool version
├── somatic_annotations
│   ├── <tumor_dna_id>.linx.breakend.tsv                # Tumor sample breakend data
│   ├── <tumor_dna_id>.linx.clusters.tsv                # Tumor sample clustered events
│   ├── <tumor_dna_id>.linx.driver.catalog.tsv          # Tumor sample driver events
│   ├── <tumor_dna_id>.linx.drivers.tsv                 #
│   ├── <tumor_dna_id>.linx.fusion.tsv                  # Tumor sample fusions
│   ├── <tumor_dna_id>.linx.links.tsv                   #
│   ├── <tumor_dna_id>.linx.neoepitope.tsv              #
│   ├── <tumor_dna_id>.linx.svs.tsv                     #
│   ├── <tumor_dna_id>.linx.vis_copy_number.tsv         #
│   ├── <tumor_dna_id>.linx.vis_fusion.tsv              #
│   ├── <tumor_dna_id>.linx.vis_gene_exon.tsv           #
│   ├── <tumor_dna_id>.linx.vis_protein_domain.tsv      #
│   ├── <tumor_dna_id>.linx.vis_segments.tsv            #
│   ├── <tumor_dna_id>.linx.vis_sv_data.tsv             #
│   └── linx.version
└── somatic_plots
    ├── all
    │   └── <tumor_dna_id>.*.png # All cluster plots
    └── reportable
        └── <tumor_dna_id>.*.png # Driver cluster plots
```

### RNA transcript analysis

**ISOFOX**

```shell
<group_id>/isofox/
├── <tumor_rna_id>.isf.alt_splice_junc.csv # Alternative splice junctions
├── <tumor_rna_id>.isf.fusions.csv         # Fusions, unfiltered
├── <tumor_rna_id>.isf.gene_collection.csv # Gene-collection fragment counts
├── <tumor_rna_id>.isf.gene_data.csv       # Gene fragment counts
├── <tumor_rna_id>.isf.pass_fusions.csv    # Fusions, filtered
├── <tumor_rna_id>.isf.retained_intron.csv # Retained introns
├── <tumor_rna_id>.isf.summary.csv         # Analysis summary
└── <tumor_rna_id>.isf.transcript_data.csv # Transcript fragment counts
```


### Oncoviral detection

**VIRUSBreakend**: Viral content and integration calling

```shell
<group_id>/virusbreakend/
├── <tumor_dna_id>.virusbreakend.vcf             # VCF with viral integration sites
└── <tumor_dna_id>.virusbreakend.vcf.summary.tsv # Analysis summary
```

**VirusInterpreter**: Post-processing

```shell
<group_id>/virusinterpreter/
└── <tumor_dna_id>.virus.annotated.tsv # Processed oncoviral call/annotation data
```

### Immune analysis

**LILAC**: HLA typing

```shell
<group_id>/lilac/
├── <tumor_dna_id>.lilac.candidates.coverage.tsv # Coverage of high scoring candidates
├── <tumor_dna_id>.lilac.qc.tsv                  # QC file
└── <tumor_dna_id>.lilac.tsv                     # Analysis summary
```

**NEO**: Neo-epitope prediction

```shell
<group_id>/neo/
├── <tumor_dna_id>.lilac.candidates.coverage.tsv # Coverage of high scoring candidates
├── <tumor_dna_id>.lilac.qc.tsv                  # QC file
└── <tumor_dna_id>.lilac.tsv                     # Analysis summary
```

### Mutational signature fitting

**SIGS**

```shell
sigs/
├── <tumor_dna_id>.sig.allocation.tsv
└── <tumor_dna_id>.sig.snv_counts.csv
```

### HRD prediction

**CHORD**

```shell
<group_id>/chord/
├── <tumor_dna_id>.chord.mutation_contexts.tsv # Counts of mutation types
└── <tumor_dna_id>.chord.prediction.tsv        # HRD predictions
```

### Tissue of origin prediction

**CUPPA**

```shell
<group_id>/cuppa/
├── <tumor_dna_id>.cuppa.pred_summ.tsv # Prediction summary               
├── <tumor_dna_id>.cuppa.vis.png       # Prediction visualisation         
├── <tumor_dna_id>.cuppa.vis_data.tsv  # Prediction visualisation raw data
└── <tumor_dna_id>.cuppa_data.tsv.gz   # Input features                   
```

### Summary report

**ORANGE**

```shell
<group_id>/orange/
├── <tumor_dna_id>.orange.pdf # Results of all tools as a PDF
└── <tumor_dna_id>.orange.json # Result raw data
```

## Work directory

When running Oncoanalyser, a work directory (default: `<current_directory>/work/`) is created that contains the input files, output 
files, and run logs for a particular tool. Once the tool is done running, the output files are 'published' (copied) to the final output 
directory. 

The work directory has the below structure:

```shell
work/
├── 06
│   └── e6f7613f50bdca27662f3d256c09e1
├── 0a
│   └── 9acb05051afef00264593f36058180
├── 1a
│   └── 9997df2e2e9978ec24b5f8e8a7bb3c
...
```

The subdirectory names are hashes and correspond to those shown in the console when running Oncoanalyser. For example, `0a/9acb05` as shown 
below is shorthand for `work/06/e6f7613f50bdca27662f3d256c09e1` as shown above, and corresponds to the `COBALT_PROFILING:COBALT` process.

> [!TIP]
> Use Tab to auto-complete directory names when navigating the work directory

```shell
...
executor >  local (28)
[-        ] process > NFCORE_ONCOANALYSER:WGTS:READ_ALIGNMENT_DNA:FASTP                            -
[-        ] process > NFCORE_ONCOANALYSER:WGTS:READ_ALIGNMENT_DNA:BWAMEM2_ALIGN                    -
[-        ] process > NFCORE_ONCOANALYSER:WGTS:READ_ALIGNMENT_RNA:STAR_ALIGN                       -
[-        ] process > NFCORE_ONCOANALYSER:WGTS:READ_ALIGNMENT_RNA:SAMTOOLS_SORT                    -
[-        ] process > NFCORE_ONCOANALYSER:WGTS:READ_ALIGNMENT_RNA:SAMBAMBA_MERGE                   -
[-        ] process > NFCORE_ONCOANALYSER:WGTS:READ_ALIGNMENT_RNA:GATK4_MARKDUPLICATES             -
[48/aa4d5c] process > NFCORE_ONCOANALYSER:WGTS:REDUX_PROCESSING:REDUX (<group_id>_<sample_id>)     [100%] 2 of 2 ✔
[2c/2acf23] process > NFCORE_ONCOANALYSER:WGTS:ISOFOX_QUANTIFICATION:ISOFOX (<group_id>)           [100%] 1 of 1 ✔
[0a/9acb05] process > NFCORE_ONCOANALYSER:WGTS:AMBER_PROFILING:AMBER (<group_id>)                  [100%] 1 of 1 ✔
[06/e6f761] process > NFCORE_ONCOANALYSER:WGTS:COBALT_PROFILING:COBALT (<group_id>)                [100%] 1 of 1 ✔
[7c/828af1] process > NFCORE_ONCOANALYSER:WGTS:ESVEE_CALLING:ESVEE_PREP (<group_id>)               [100%] 1 of 1 ✔
[e1/182433] process > NFCORE_ONCOANALYSER:WGTS:ESVEE_CALLING:ESVEE_ASSEMBLE (<group_id>)           [100%] 1 of 1 ✔
[76/0da3ee] process > NFCORE_ONCOANALYSER:WGTS:ESVEE_CALLING:ESVEE_DEPTH_ANNOTATOR (<group_id>)    [100%] 1 of 1 ✔
[41/49f1f8] process > NFCORE_ONCOANALYSER:WGTS:ESVEE_CALLING:ESVEE_CALL (<group_id>)               [100%] 1 of 1 ✔
[ce/0f6b20] process > NFCORE_ONCOANALYSER:WGTS:SAGE_CALLING:GERMLINE (<group_id>)                  [100%] 1 of 1 ✔
[5e/be6aab] process > NFCORE_ONCOANALYSER:WGTS:SAGE_CALLING:SOMATIC (<group_id>)                   [100%] 1 of 1 ✔
[45/88540d] process > NFCORE_ONCOANALYSER:WGTS:PAVE_ANNOTATION:GERMLINE (<group_id>)               [100%] 1 of 1 ✔
[e2/279465] process > NFCORE_ONCOANALYSER:WGTS:PAVE_ANNOTATION:SOMATIC (<group_id>)                [100%] 1 of 1 ✔
[ff/37883b] process > NFCORE_ONCOANALYSER:WGTS:PURPLE_CALLING:PURPLE (<group_id>)                  [100%] 1 of 1 ✔
[d0/7ebc71] process > NFCORE_ONCOANALYSER:WGTS:SAGE_APPEND:GERMLINE (<group_id>)                   [100%] 1 of 1 ✔
[1c/0b3f55] process > NFCORE_ONCOANALYSER:WGTS:SAGE_APPEND:SOMATIC (<group_id>)                    [100%] 1 of 1 ✔
[87/0118e3] process > NFCORE_ONCOANALYSER:WGTS:LINX_ANNOTATION:GERMLINE (<group_id>)               [100%] 1 of 1 ✔
[1a/9997df] process > NFCORE_ONCOANALYSER:WGTS:LINX_ANNOTATION:SOMATIC (<group_id>)                [100%] 1 of 1 ✔
[a8/22db2b] process > NFCORE_ONCOANALYSER:WGTS:LINX_PLOTTING:VISUALISER (<group_id>)               [100%] 1 of 1 ✔
[dc/da6010] process > NFCORE_ONCOANALYSER:WGTS:BAMTOOLS_METRICS:BAMTOOLS (<group_id>_<sample_id>)  [100%] 2 of 2 ✔
[b5/5c54f6] process > NFCORE_ONCOANALYSER:WGTS:SIGS_FITTING:SIGS (<group_id>)                      [100%] 1 of 1 ✔
[71/701751] process > NFCORE_ONCOANALYSER:WGTS:CHORD_PREDICTION:CHORD (<group_id>)                 [100%] 1 of 1 ✔
[bc/6191b2] process > NFCORE_ONCOANALYSER:WGTS:LILAC_CALLING:LILAC (<group_id>)                    [100%] 1 of 1 ✔
[51/153ee1] process > NFCORE_ONCOANALYSER:WGTS:VIRUSBREAKEND_CALLING:VIRUSBREAKEND (<group_id>)    [100%] 1 of 1 ✔
[88/fee470] process > NFCORE_ONCOANALYSER:WGTS:VIRUSBREAKEND_CALLING:VIRUSINTERPRETER (<group_id>) [100%] 1 of 1 ✔
[28/6e9733] process > NFCORE_ONCOANALYSER:WGTS:CUPPA_PREDICTION:CUPPA (<group_id>)                 [100%] 1 of 1 ✔
[e0/2e5797] process > NFCORE_ONCOANALYSER:WGTS:ORANGE_REPORTING:ORANGE (<group_id>)                [100%] 1 of 1 ✔
...
```

Below is an example of the contents of the `COBALT_PROFILING:COBALT` process work directory. 

```shell
work/06/
└── e6f7613f50bdca27662f3d256c09e1
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .command.trace
    ├── .exitcode
    ├── <normal_dna_id>.redux.bam ->  /path/to/work/32/6d0191b876479d1a0c3c4a4c39733d/<normal_dna_id>.redux.bam
    ├── <normal_dna_id>.redux.bam.bai ->  /path/to/work/32/6d0191b876479d1a0c3c4a4c39733d/<normal_dna_id>.redux.bam.bai
    ├── <tumor_dna_id>.redux.bam ->  /path/to/work/48/aa4d5cecc431bfe3fef5e85d922272/<tumor_dna_id>.redux.bam
    ├── <tumor_dna_id>.redux.bam.bai ->  /path/to/work/48/aa4d5cecc431bfe3fef5e85d922272/<tumor_dna_id>.redux.bam.bai
    ├── GC_profile.1000bp.37.cnp -> /path/to/hmftools/dna/copy_number/GC_profile.1000bp.37.cnp
    ├── cobalt
    │   ├── <normal_dna_id>.cobalt.gc.median.tsv
    │   ├── <normal_dna_id>.cobalt.ratio.median.tsv
    │   ├── <normal_dna_id>.cobalt.ratio.pcf
    │   ├── <tumor_dna_id>.cobalt.gc.median.tsv
    │   ├── <tumor_dna_id>.cobalt.ratio.pcf
    │   ├── <tumor_dna_id>.cobalt.ratio.tsv.gz
    │   └── cobalt.version
    └── versions.yml
```

Tool work directories have a consistent structure:
- `.command.sh`: Bash command used to run the tool <ins>within the Docker/Singularity container<ins>
- `.command.log`,  `.command.err`, `.command.out`: Run logs
- `versions.yml`: Tool version
- Tool outputs generally are written to a directory of the same name (e.g. `cobalt/`)
- Input files are symlinked into the tool work directory (e.g. `<tumor_dna_id>.redux.bam ->  ...`). This is done so that under the hood the 
tool work directory can simply be [mounted](https://docs.docker.com/engine/storage/bind-mounts/) within the container.

## Acknowledgements

Oncoanalyser was written by Stephen Watts at the [University of Melbourne Centre for Cancer
Research](https://mdhs.unimelb.edu.au/centre-for-cancer-research) with the support of Oliver Hofmann and the Hartwig
Medical Foundation Australia.
