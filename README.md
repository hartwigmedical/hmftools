# HMF Tools

This repository contains the suite of tools used in the Hartwig Medical Foundation whole genome,targeted DNA and whole transcriptome analysis pipeline.  

## DNA Pipeline Components

![HMF_Pipeline](./pipeline/hmf_tools_pipeline.png)

The table below has links for the full functional detail and configuration for the key components in our DNA pipeline. 

The versions match those used in the current HMF GCP pipeline which can be run using [Platinum](https://github.com/hartwigmedical/platinum). 

Demo scripts to run each component in turn are provided [here](./pipeline/scripts/).

Component | Description | Current Version
---|---|---
[Sage](./sage/README.md) | Point mutation variant calling and filtering | [3.0](https://github.com/hartwigmedical/hmftools/releases/tag/sage-v3.0)
[Pave](./pave/README.md) | Point mutation annotation and gene impact | [1.2](https://github.com/hartwigmedical/hmftools/releases/tag/pave-v1.2)
[Amber](./amber/README.md) | Generate a tumor BAF file for Purple's copy number fit | [3.9](https://github.com/hartwigmedical/hmftools/releases/tag/amber-v3.9)
[Cobalt](./cobalt/README.md) | Determines the read depth ratios for Purple's copy number fit | [1.13](https://github.com/hartwigmedical/hmftools/releases/tag/cobalt-v1.13)
[Gripss](./gripss/README.md) | SV filtering | [2.1](https://github.com/hartwigmedical/hmftools/releases/tag/gripss-v2.1) 
[Purple](./purple/README.md) | Estimates copy number, purity and ploidy, and identifies driver events | [3.4](https://github.com/hartwigmedical/hmftools/releases/tag/purple-v3.4)
[Linx](./linx/README.md) | SV annotation, clustering & chaining, fusion and disruption calling | [1.19](https://github.com/hartwigmedical/hmftools/releases/tag/linx-v1.19)
[Lilac](./lilac/README.md) | HLA typing | [1.2](https://github.com/hartwigmedical/hmftools/releases/tag/lilac-v1.2)
[Teal](./lilac/README.md) | Measures telomere content and estimates telomeric length | [1.0.1](https://github.com/hartwigmedical/hmftools/releases/tag/teal-v1.0.1)
[Cuppa](./cuppa/README.md) | Tissue of origin prediction from WGS/WTS | [1.6](https://github.com/hartwigmedical/hmftools/releases/tag/cuppa-v1.6)

Note that [Chord](https://github.com/UMCUGenetics/CHORD) and [GRIDSS](https://github.com/PapenfussLab/gridss) are supported externally to HMF.

## Actionability and Clinical Reporting Tools

Component | Description | Current Version
---|---|---
[Serve](./serve/README.md) | Harmonisation of evidence from clinical annotation databases | [1.11](https://github.com/hartwigmedical/hmftools/releases/tag/serve-v1.11) 
[Protect](./protect/README.md) | Matching of molecular results to treatments and clinical trials | [2.2](https://github.com/hartwigmedical/hmftools/releases/tag/protect-v2.2)
[Rose](./rose/README.md) | Actionability of clinically relevant molecular findings | [1.0](https://github.com/hartwigmedical/hmftools/releases/tag/rose-v1.0)
[Virus Interpreter](./virus-interpreter/README.md) | Filtering, annotation and interpretation of virus breakend data | [1.2](https://github.com/hartwigmedical/hmftools/releases/tag/virus-interpreter-v1.2)
[Orange](./orange/README.md) | PDF summary report and JSON file of all WGS output | [1.9](https://github.com/hartwigmedical/hmftools/releases/tag/orange-v1.9)

## RNA Tools

Component | Description | Current Version
---|---|---
[Isofox](./isofox/README.md) | WTS Transcript Abundance, Fusions & Novel Splice Junctions | [1.5](https://github.com/hartwigmedical/hmftools/releases/tag/isofox-v1.5)

## Resource Files
Resources used by the HMF algorithms can be downloaded from https://resources.hartwigmedicalfoundation.nl

The complete set of resource files required to run the DNA pipeline described above can be found here.

