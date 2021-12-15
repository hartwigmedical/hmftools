# GRIDSS Post Somatic Software (GRIPSS)

This version has been superceded by:
https://github.com/hartwigmedical/hmftools/blob/master/gripss/README.md


## Usage for version 1.12 and earlier

```
java -Xms4G -Xmx16G -cp gripss.jar com.hartwig.hmftools.gripsskt.GripssApplicationKt \
   -tumor SAMPLE_T \
   -reference SAMPLE_N \
   -ref_genome /path/to/Homo_sapiens_assembly.fasta \
   -breakend_pon /path/to/gridss_pon_single_breakend.bed \
   -breakpoint_pon /path/to/gridss_pon_breakpoint.bedpe \
   -breakpoint_hotspot /path/to/KnownFusionPairs.bedpe \
   -input_vcf /path/to/SAMPLE_T.gridss.unfiltered.vcf.gz \
   -output_vcf /path/to/SAMPLE_T.gripss.vcf.gz 
```

Run a second time to remove all SVs except PASS and PON:
```
java -Xms4G -Xmx16G -cp gripss.jar com.hartwig.hmftools.gripsskt.GripssHardFilterApplicationKt \
   -input_vcf /path/to/SAMPLE_T.gripss.vcf.gz \
   -output_vcf /path/to/SAMPLE_T.gripss.filtered.vcf.gz 
```

These two files are used in purple as the structural variant recovery vcf and structural variant vcf respectively.


## Version History and Download Links
- [1.11](https://github.com/hartwigmedical/hmftools/releases/tag/gripss-v1.11)
  - Fix hotspot matching
  - Replace other ambiguous bases in bam with N
  - Fixed description of TAF: Tumor allelic frequency (fragment support / total support)
  - Do not attempt transitive linking if there is 500,000+ variants
  - Changed default value of hardMaxNormalRelativeSupport from 0.06 to 0.08 
- [1.10](https://github.com/hartwigmedical/hmftools/releases/tag/gripss-v1.10)
  - Removed viral insertion exception from normal relative support filter
  - Add EventType flag [DEL, INS, DUP, INV, SGL, BND]
  - Check IC flag whenever we check SR 
- [1.9](https://github.com/hartwigmedical/hmftools/releases/tag/gripss-v1.9)
  - Warn if malformed BEALN field in input
  - Only include specified samples in VCF output
  - `tumor` argument is now mandatory
  - `reference` argument is still optional but excluding it will run in [tumor-only mode](#tumor-only-mode)
  - HardMaxNormalAbsoluteSupport does not trigger unless SoftMaxNormalRelativeSupport does also
- [1.8](https://github.com/hartwigmedical/hmftools/releases/tag/gripss-v1.8)
  - Legs of variant always have same filters
- [1.7](https://github.com/hartwigmedical/hmftools/releases/tag/gripss-v1.7)
  - Fix bug in BEALN interpretation
- [1.6](https://github.com/hartwigmedical/hmftools/releases/tag/gripss-v1.6)
  - Fix bug where alt does not always realign
- [1.5](https://github.com/hartwigmedical/hmftools/releases/tag/gripss-v1.5)
  - Tumor only support
  - Added GripssHardFilterApplicationKt application
- [1.4](https://github.com/hartwigmedical/hmftools/releases/tag/gripss-v1.4)
  - Fix bug when trying to use non-default parameter values
- [1.3](https://github.com/hartwigmedical/hmftools/releases/tag/gripss-v1.3)
  - Fixed bug to allow output VCF name to be relative path
  - Fixed bug handling contigs with colons eg HLA-A*23:01:01
- [1.2](https://github.com/hartwigmedical/hmftools/releases/tag/gripss-v1.2)
  - Hotspots cannot recover HardMinQualFiltered variants
  - Ignore maxNormalRelativeSupport filter for single breakends with viral sequence alignments
  - BEALN, INSRMRC and INSRMRT are now mandatory requirements
- [1.1](https://github.com/hartwigmedical/hmftools/releases/tag/gripss-v1.1)
  - Added HardMinTumorQual [100] filter
  - Added HardMaxNormalAbsoluteSupport [3] filter 
  - Added HardMaxNormalRelativeSupport [0.06] filter to replace MaxNormalRelativeSupport [0.03]
  - Added SoftMaxNormalRelativeSupport [0.03] filter
  - Qual filters and QUAL field in output VCF now refer to tumor qual only
  - Added optional `reference` and `tumor` parameters
  - Rescue DSB linked mobile element insertions
- [1.0](https://github.com/hartwigmedical/hmftools/releases/tag/gripss-v1.0)
  - Initial Release 
