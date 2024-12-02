# GRIDSS Post Somatic Software (GRIPSS)

GRIPSS applies a set of filtering and post processing steps on GRIDSS paired tumor-normal output to produce a high confidence set of somatic SV for a tumor sample.
GRIPSS processes the GRIDSS output and produces a somatic vcf.

## Tumor-only / Germline mode
The `reference` argument is optional and if not supplied, GRIPSS will run in 'tumor-only' mode in which case  all filters that require the normal sample are de-activated. 
This includes `minNormalCoverage`, `minRelativeCoverage`, `maxNormalSupport`, `shortSRNormalSupport`.
Single breakends are not called in tumor only mode as there are many germline artefacts. GRIPSS can be run in a germline mode by setting the 'sample' to the germline and not supplying a 'reference' argument.
In the hartwig pipeline we use the same filter settings for tumor and germline structural variants in GRIDSS (see below). 
 
# Algorithm

There are 4 steps in GRIPSS described in detail below:
  1. [Hard filters](#1-hard-filters)
  2. [Realignment](#2-realignment)
  3. [Soft filters](#3-soft-filters)
  4. [Linkage, deduplication and rescue](#4-linkage-deduplication-and-rescue)

## 1. Hard filters

Three hard filters are applied upfront before other processing occurs:
* NO_MATE - Any non single breakend with no mate is filtered
* MINIMUM_TUMOR_QUAL - Any variant with QUAL < 100 is filtered
* MAX_NORMAL_SUPPORT - Any variant with normalSupport > min(max(3, 3% * tumorSupport), 8% * tumorSupport) is filtered as likely germline or artefact unless it links a pair of genes in the known pathogenic fusion list via translocation or local break junction of length more than 10kb. Ideally we would not allow any support for the variant in the normal, but contamination of the blood with tumor DNA is not uncommon.

Type | Config | Default 
---|---|---
NO_MATE | N/A | N/A  
MINIMUM_TUMOR_QUAL | hard_min_tumor_qual | 100  
MAX_NORMAL_SUPPORT | hard_max_normal_absolute_support | 3
MAX_NORMAL_SUPPORT | hard_max_normal_relative_support | 0.08
MAX_NORMAL_SUPPORT | soft_max_normal_relative_support* | 0.03

<nowiki>*</nowiki> This filter is also used as a soft filter

## 2. Realignment

We realign imprecise variants where the breakend can be precisely resolved but the length of the insert sequence is unknown.  By default GRIDSS offsets these variants by the uncertainty of the insert sequence with a wide CIPOS.  GRIPSS realigns the variant to the earliest possible base in the uncertainty window which is the most likely base for the soft clipping.

For the purposes of backwards compatibility we also perform 3 other fixes to correct errors in earlier GRIDSS versions
* Breakends are shifted to the centre of homology.  inexact homology bounds are also realigned
* Ensure that breakends pairs are internally consistent in their positions
* Ensure that local and remote inexact homology are internally consistent.

## 3. Soft filters
 
The following filters are applied to all variants that aren't hotspots:

Filter | Config | Default | Description / purpose
---|---|---|---
minQual | min_qual_break_point, min_qual_break_end | 400, SGL = 500 | Minimum absolute tumor support for variant. For PMS2 (an important MMR gene) only (including 10kb upstream and downstream) the minimum qual is set to half the normal value due to the poor mappability of PMS2 due to it's close homolog PMS2CL.
minNormalCoverage | min_normal_coverage | 8 | Variants with low coverage in germline may be germline variants.
minTumorAF | min_tumor_af | 0.005 | Low AF variants in high depth regions may be artefacts. SGLs use S0.015.
imprecise | N/A | FALSE | Imprecise variants may be artefacts linking low mappability regions of the genome.   
PON | PON files | FALSE | Breakpoint must be found < 3 times in our cohort in pon_sv_file (or pon_sgl_file in the case of breakends). The Hartwig PON files are generated from ~3800 germline samples. The PON excludes imprecise calls and breakpoints <75 qual score and breakends < 428 qual score.  Inexact homology is allowed in overlap and an additional 3p margin of error is allowed for. 
maxPolyAHomLength | N/A | 6 | Variants with long poly-A homology are frequent artefacts at low VAF
maxPolyGLength | N/A | 16 | Long stretches of poly-G/poly-C are extremely rare in the ref genome but are known sequencer artefacts.  Single breakends with insert sequences containing long polyG homopolymers are filtered.   This filter is also applied to break junctions where 1 end maps in any of the following POLY-G regions (v38: {chr2:32,916,190-32,916,630; chr4:41,216,410-41,216,450; chr17:44,569,050-44,569,090}; v37: {2:33,141,260-33,141,700; 4:41,218,427-41,218,467; 17:42646418-42646458}).
maxNormalRelativeSupport | soft_max_normal_relative_support* | 0.03 | Too many support reads from the normal sample relative to the tumor. Indicates the variant is likely germline or artefact.
qualPerAD | qual_per_ad | 30 | Require an average qual contributrion per VF support.  ly in targeted mode.  Does not apply to HOTSPOTS

<nowiki>*</nowiki> This argument is also used for one of the hard filters

We also have 9 special filters applying to specific variant categories:

Filter | Config | Default | Scope | Description 
---|---|---|---|---
minLength | min_length | 32 | DEL, DUP & INS | Minimum absolute length (including insert sequence length) for short DEL and DUP SV to be called. 
minSingleInsertLength | N/A | 16 | SGL | Minimum insert sequence length for a single breakend
singleStrandBias | max_short_strand_bias | 0.05<SB<0.95 | SGL (excluding polyA tails) | Minimum/maximum proportion of reads from the forward strand supporting the single breakend.   PMS2 is excluded from this rule due to the poor mappability and high homology
maxHomLengthShortInv | max_hom_length_short_inv | 6 | INV(<50b) | Very short INV with high homology are a common sequencer artefact
discordantPairSupport | N/A | TRUE | INV(<50b) | Breakpoints must have at least 1 read mapped at each end.  Filters if ASRP=RP=0.
shortStrandBias | N/A | TRUE | INS,DEL & DUP(<1kb) | Short DEL and DUP must be strand balanced
shortSRTumorSupport | N/A | TRUE | INS,DEL & DUP(<1kb) | Short DELs and DUPs must be supported by at least 1 split read or in the case of very short DEL and INS at least 1 supporting indel containing read.
shortSRNormalSupport | N/A | FALSE | INS,DEL & DUP(<1kb) | Short DELs and DUPs must not be supported by 1 split read or 1 indel containing read in the normal 
shortDelInsArtefact | N/A | TRUE | DEL(<1kb) | Filter any short DEL where the insert sequence length + 1 = deletion length, unless the insert sequence is identical to the reverse complement (ie a short reciprocal inversion).  This is a known GRIDSS artefact.

For targeted panel data filter are customised and single breakends are not called. Please see [here](https://github.com/hartwigmedical/hmftools/blob/master/pipeline/README_TARGETED.md).

## 4. Linkage, deduplication and rescue

### A. Assembly linkage

Variants located on a single assembly are given a unique 'asm' identifier.

### B. Deduplication and transitive linkage

The GRIDSS output may contain structural variants which may be duplicated either by a single SV or by a chain of SVs with breakends proximate to each other which may or may not already be linked by assembly.
In the case where a variant is duplicated by a chain, we term this variant the spanning variant and these links to be transitive links.

For a variant to be marked as a duplicate, we must find 2 candidate transitive breakends which match the orientation and position of the spanning variant, within CIPOS bounds and allowing for the insert sequence length.
The CIPOS bounds for imprecise variants are ignored for this purpose and the match must be exact (see realignment section above).    

The candidate transitive breakends must be linkable in a continuous chain as one of the following cases:

* same variant - opposite breakends of the same SV
* same assembly - the 2 transitive SVs are part of the same assembly and oriented away from each other
* 1 transitive jump - the far breakend of the 2 transitive SVs / assemblies face each other and each link must be less than 1000 bases
* 2 transitive jumps - the far breakend of the 2 transitive SVs / assemblies both face opposite ends of a 3rd SV or assembly and each link must be less than 1000 bases

As a spanning variant may have multiple alternate paths, we first consider only assembly linked paths favouring those with the fewest jumps.  
If an assembly-only solution is unavailable we include up to two transitive jumps, but this must result in a single alternate path. 
If there are multiple alternate paths with transitive links, none will be selected.
 
If the deduplicated spanning variant is PRECISE, then the length of the insert sequence of the spanning variant must match the entire chain length of the transitive variants (again allowing for CIPOS bounds and insert sequence length of precise variants).  

Any single breakend which matches the position and orientation of another breakend or breakjunction (within CIPOS bounds) is also filtered as DEDUP.
GRIPSS prioritises retaining breakends / breakjunctions linked by assembly or transitive link, then breakends that are passing and finally by highest qual score. 

### C. Linkage by double stranded break

Double stranded break sites can lead to 2 proximate breakends in very close proximity with opposite orientation and either a small gap in between or a small overlap, frequently <30 bases.   GRIPSS links breakends with a unique 'DSB' id  where there is one and only one (non-DEDUPed) breakend within 30 bases (allowing for CIPOS bounds) with opposite orientation and both breakends have qual > 100.

### D. Rescue

Any breakend that is linked to a PASS breakend (by one of the 3 above rules) and is NOT filtered as DEDUP is rescued from soft filtering and marked as PASS. 
Breakend pairs that link a pair of genes to make a known pathogenic fusions are also rescued for all soft filters except maxPolyAHomLength.

To improve detection of mobile element insertions, we also rescue pairs of breakends or breakjunctions which are linked by ‘DSB’ and NOT PON filtered, 
with combined qual > 500 and with at least one of the breakends having the characteristic poly-A insert sequence tail of a mobile element insertion. 
We define a poly-A tail as 16 of the last 18 bases of the insert sequence are A. At the insertion site, negative oriented breakends must have poly-A tails 
at the end of the insert sequence and positive oriented breakends must have poly-T at the start of the insert sequence (if inserted on the reverse strand).

Note that for DSB and hotspot rescue, neither the rescued variant nor the rescuing variant is permitted to be a DEL, INS or DUP < 10kb in length.  

### E. Repeat masker annotation

Each variant with a BEALN (insert sequence alignment) is annotated with 4 additional fields from the entry with the most overlapping sequence from the repeat masker hg19.fa.out.gz file  which covers at least 10% and 20 bases of the alignment:

* INSRMP= Portion of inserted sequence whose alignment overlaps the repeatmasker repeat. 1.0 indicates the inserted sequence entirely mapping to the repeat
* INSRMRC=Inserted sequence repeatmasker repeat class  [matching repeat]
* INSRMRO=Inserted sequence repeatmasker repeat orientation
* INSRMRT=Inserted sequence repeatmasker repeat type [repeat class/family]

Note in the special case where more than 90% of the candidate alignment bases are  A or T then this is just reported as a simple_repeat A(n) or T(n) regardless of the repeatmasker annotation.

The INSRMRC and INSRMRT are used in various clustering rules downstream in LINX.   This functionality replaces fully the previous repeat masker post process for GRIDSS.  


## Counting Conventions in GRIPSS
- **Fragment support** - the count of supporting fragments for breakpoints is set to VF field from GRIDSS and for single breakends is set to BVF (with an exception that if a single breakend has BSC=BASRP=BASSR=0 then read support is set to 0).   The fragment support is used in the various support filters as described above
- **Qual score** - for breakpoints the qual is set to the qual.   For breakends the BAQ field (ie sum of qual for assembled reads supporting the breakend) is used for the qual except where the insert sequence has a poly-A tail in which case the BQ (qual of all reads supporting the breakend) is used.

