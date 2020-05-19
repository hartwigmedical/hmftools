# GRIDSS Post Somatic Script (GRIPSS)

GRIPSS applies a set of filtering and post processing steps on GRIDSS paired tumor-normal output to produce a high confidence set of somatic SV for a tumor sample.    GRIPSS inputs the raw GRIDSS vcf and outputs a somatic vcf.


# Algorithm

There are 6 key steps in GRIPSS described in detail below:
  1. [Hard Filters](#1-hard-filters)
  2. [Realignment](#2-realignment)
  3. [Soft Filters](#3-soft-filters)
  4. [Linkage, Deduplication and Rescue](#4-linkage,-deduplication-and-rescue)
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
 
<TO DO>

## 4. Linkage, deduplication and rescue

### A. Assembly linkage

<TO DO>

### B. Transitive linkage and deduplication

<TO DO>

### C. Linkage by double stranded break

<TO DO>

### D. Rescue

Any breakend that is linked to a PASS breakend (by one of the 3 above rules) is rescued from soft filtering and marked as PASS.    Breakend pairs that link a pair of genes to make a known pathogenic fusions are also rescued regardless of soft filtering.

<TO DO - specify known file>

## 5. PON Filtering

<TO DO>

# Version History and Download Links
 - Upcoming
