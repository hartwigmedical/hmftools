# HMF Tools
This repository contains the suite of tools used in the Hartwig Medical Foundation whole genome and transcriptome analysis pipeline.

![HMF_Pipeline](pipeline/hmf_tools_pipeline.png)

## DNA Pipeline Components
Click each component's link for the full functional detail and configuration.

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

### Resource Files
Resources used by the HMF algorithms can be downloaded from https://resources.hartwigmedicalfoundation.nl

The complete set of resource files required to run the DNA pipeline described above can be found here.

## Further Tools
- [Isofox](./isofox/README.md)
- [Serve](./serve/README.md)
- [Protect](./protect/README.md)
- [Cuppa](./cuppa/README.md)
- [Virus Interpreter](./virus-interpreter/README.md)
- [Rose](./rose/README.md)
- [Orange](./orange/README.md)
