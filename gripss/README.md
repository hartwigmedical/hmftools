# GRIDSS Post Somatic Script (GRIPSS)

GRIPSS applies a set of filtering and post processing steps on GRIDSS paired tumor-normal output to produce a high confidence set of somatic SV for a tumor sample.    GRIPSS inputs the raw GRIDSS vcf and outputs a somatic vcf.


# Algorithm

There are 6 key steps in GRIPSS described in detail below:
  1. [Hard Filters](#1-hard-filters)
  2. [Realignment](#2-realignment)
  3. [Soft Filters](#3-soft-filters)
  4. [Linkage, Deduplication and Rescue](#4-linkage-deduplication-and-rescue)
  5. [Pon Filtering](#5-pon-filtering)

## 1. Hard filters

Two hard filters are applied upfront before other processing occurs:
* NO_MATE - Any non single breakend with no mate is filtered
* MAX_NORMAL_SUPPORT - Any variant with normalSupport > 3% of tumor support is filtered as likely germline or artefact.  Ideally we would not allow any support for the variant in the normal, but contamination of the blood with tumor DNA is not uncommon.

## 2. Realignment

We realign imprecise variants where the breakend can be precisely resolved but the length of the insert sequence is unknown.  By default GRIDSS offsets these variants by the uncertainty of the insert sequence with a wide CIPOS.

For the purposes of backwards compatibility we also perform 3 other fixes to correct errors in earlier GRIDSS versions
* Breakends are shifted to the centre of homology.  
* Ensure that breakends pairs are internally consistent in their positions
* Ensure that local and remote inexact homology are internally consistent.

## 3. Soft filters
 
The following filters are applied to all variants

Filter | Default | Description / purpose
---|---|---
minQual | breakpoints: 350; single breakend:1000 | Minimum absolute tumor support for variant
minNormalCoverage | 8 | Variants with low coverage in germline may be germline variants.
minTumorAF | 0.5 | Low AF variants in high depth regions may be artefacts
imprecise | FALSE | Imprecise variants may be artefacts linking low mappability regions of the genome  
maxInexactHomologyLength | 50 | Very long inexact homology may also be artefacts linking low mappability regions of the genome
DPsupport | TRUE | Variants (except for DEL and DUP < 1000 bases) must have at least 1 read mapped at each end.   Avoids artefacts linking regions of low mapability.   Not suitable for non paired reads or very short fragment sizes.  

Single breakends have 2 additional filters:

Filter | Default | Description 
---|---|---
breakendAssemblyReadPair | TRUE | ???
maxPolyGLength | 16 | Long stretches of polyG/polyC are extremely rare in the ref genome but are known sequencer artefacts

We also have 7 special filters applying to specific short variant categories:

Filter | Default | Scope | Description 
---|---|---|---
minDelDupLength | 32 | DEL, DUP & INS | Minimum absolute length (including insert sequence length) for short DEL and DUP SV to be called. 
shortInvMaxHomLength | 6 | INV(<40b) | 
shortDelMaxHomLength | 6 | DEL(<1kb) | 
shortDelDupStrandBias | true | DEL & DUP(<1kb) | Enable base quality recalibration
shortDelDupSRSupport | TRUE | DEL & DUP(<1kb) | Short DELs and DUPs must be supported by at least 1 split read
shortDelDupNormalSRSupport | FALSE | DEL & DUP(<1kb) | Short DELs and DUPs must not be supported by 1 split read in the normal
small.replacement.fp | | DEL(<1kb) | *TO DO!!!*




## 4. Linkage, deduplication and rescue

### A. Assembly linkage

<Jon - can you write

### B. Deduplication and transitive linkage

The GRIDSS output may contain structural variants which may be duplicated either by a single SV or by a chain of SVs with breakends proximate to each other which may or may not already be linked by assembly.   In the case where a variant is duplicated by a chain, we term this variant the spanning variant and these links to be transitive links.

For a variant to be marked as a duplicate, we must find 2 candidate transitive breakends which match the orientation and position (within CIPOS bounds) of the spanning variant. The candidate transitive breakends must be linkable in a continuous chain as one of the following cases:

* same variant - opposite breakends of the same SV
* same assembly - the 2 transitive SVs are part of the same assembly and oriented away from each other
* 1 transitive jump - the far breakend of the 2 transitive SVs / assemblies face each other and each link must be less than 1000 bases
* 2 transitive jumps - the far breakend of the 2 transitive SVs / assemblies both face opposite ends of a 3rd SV or assembly and each link must be less than 1000 bases

Further, all the transitive links must be precise variants.   If the deduplicated spanning variant is PRECISE, then the length of the insert sequence of the spanning variant must match the entire chain length of the transitive variants (within CIPOS bounds).  

Any single breakend which matches the position and orientation of another breakend (within CIPOS bounds) is also filtered as DEDUP (if the matching variant is also a single breakend then the highest scoring passing breakend will be kept).

### C. Linkage by double stranded break

<TO DO>

### D. Rescue

Any breakend that is linked to a PASS breakend (by one of the 3 above rules) is rescued from soft filtering and marked as PASS.    Breakend pairs that link a pair of genes to make a known pathogenic fusions are also rescued regardless of soft filtering.

<TO DO - specify known file>

## 5. PON Filtering

<TO DO>
  
