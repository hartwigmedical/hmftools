# HMF Tools

This repository contains the suite of tools used in the Hartwig Medical Foundation whole genome,targeted DNA and whole transcriptome analysis pipeline.  

## DNA Pipeline Components

![HMF_Pipeline](./pipeline/hmf_tools_pipeline.png)

The table below has links for the full functional detail and configuration for each component. The versions match those used in the current HMF GCP pipeline which can be run using [Platinum](https://github.com/hartwigmedical/platinum).

#### Current versions
Component | Description                                                            | Current Version
---|------------------------------------------------------------------------|---
[Amber](./amber/README.md) | Generate a tumor BAF file for Purple's copy number fit                 | [3.9](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v3.9)
[Cobalt](./cobalt/README.md) | Determines the read depth ratios for Purple's copy number fit          | [1.13](https://github.com/hartwigmedical/hmftools/releases/tag/cobalt-v1.13)
[Cuppa](./cuppa/README.md) | Tissue of origin prediction from WGS/WTS                               | [1.7](https://github.com/hartwigmedical/hmftools/releases/tag/cuppa-v1.7)
[Gripss](./gripss/README.md) | SV filtering                                                           | [2.3](https://github.com/hartwigmedical/hmftools/releases/tag/gripss-v2.3.2) 
[Lilac](./lilac/README.md) | HLA typing                                                             | [1.4](https://github.com/hartwigmedical/hmftools/releases/tag/lilac-v1.4.1)
[Linx](./linx/README.md) | SV annotation, clustering & chaining, fusion and disruption calling    | [1.22](https://github.com/hartwigmedical/hmftools/releases/tag/linx-v1.22)
[Pave](./pave/README.md) | Point mutation annotation and gene impact                              | [1.4](https://github.com/hartwigmedical/hmftools/releases/tag/pave-v1.4)
[Purple](./purple/README.md) | Estimates copy number, purity and ploidy, and identifies driver events | [3.7](https://github.com/hartwigmedical/hmftools/releases/tag/purple-v3.7.2)
[Sage](./sage/README.md) | Point mutation variant calling and filtering                           | [3.2.2](https://github.com/hartwigmedical/hmftools/releases/tag/sage-v3.2.2)
[Teal](./teal/README.md) | Measures telomere content and estimates telomeric length               | [1.0.1](https://github.com/hartwigmedical/hmftools/releases/tag/teal-v1.0.1)
[SvPrep](./sv-prep/README.md) | Pre-GRIDSS BAM filtering | [1.0](https://github.com/hartwigmedical/hmftools/releases/tag/sv-prep-v1.0.1)

The following external tools are used in the pipeline:

Component | Description | Current Version
---|---|---
[GRIDSS](https://github.com/PapenfussLab/gridss) | Structural variant calling | [2.13.2](https://github.com/PapenfussLab/gridss/releases/tag/v2.13.2)
[Chord](https://github.com/UMCUGenetics/CHORD) | Homologous Recombination Deficiency detection | [2.0](https://github.com/UMCUGenetics/CHORD/releases/tag/2.00)


### Demo pipeline
An example WGS pipeline which runs each of these components in turn is detailed [here](./pipeline/README_WGS.md).
An example targeted-panel pipeline, including support for the HMF and TSO-500 panels, is detailed [here](./pipeline/README_PANEL.md).

## Resource files
Resource files for each component (GRCh37 and GRCh38) are available to download from [HMFTools-Resources > DNA-Resources](https://console.cloud.google.com/storage/browser/hmf-public/HMFtools-Resources/dna_pipeline/). 

## Actionability and Clinical Reporting Tools

Component | Description                                                          | Current Version
---|----------------------------------------------------------------------|--- 
[Protect](https://github.com/hartwigmedical/oncoact/tree/master/protect/README.md) | Matching of molecular results to treatments and clinical trials      | [2.3](https://github.com/hartwigmedical/hmftools/releases/tag/protect-v2.3)
[Rose](https://github.com/hartwigmedical/oncoact/tree/master/rose/README.md) | Actionability of clinically relevant molecular findings              | [1.3](https://github.com/hartwigmedical/hmftools/releases/tag/rose-v1.3)
[Virus Interpreter](./virus-interpreter/README.md) | Filtering, annotation and interpretation of virus breakend data      | [1.2](https://github.com/hartwigmedical/hmftools/releases/tag/virus-interpreter-v1.2)
[Orange](./orange/README.md) | PDF summary report and JSON file of all WGS output                   | [1.10](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v1.10)
Patient-reporter | PDF summary report and JSON file of all clinical relevant WGS output | 7.25.1

## RNA Tools

Component | Description | Current Version
---|---|---
[Isofox](./isofox/README.md) | WTS Transcript Abundance, Fusions & Novel Splice Junctions | [1.5](https://github.com/hartwigmedical/hmftools/releases/tag/isofox-v1.5)

