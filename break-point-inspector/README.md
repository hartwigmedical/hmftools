# Break Point Inspector (BPI)

## Introduction

Manta by Illumina works well a majority of the time in a somatic variant detection mode.
However, based on our experience it has a significant false positive rate.

We have benchmarked against a COLO829 (50x depth) and found that a default parameter Manta run results in a 37% false positive rate (38 / 103 calls).

This has been manually verified by inspection in IGV, as well as comparison to a default Gridss run and the Conserting study. We have generated a non-exhaustive truth set for COLO829 from the results.

You can view our truth set here:
[Hartwig Medical Foundation - COLO829 Truth Set](https://docs.google.com/spreadsheets/d/e/2PACX-1vTF5IeIoQXz-Dny0eauDbTDtyIi2nL8fTKcLr_ByHO2BOClxrj3SQ-GJBRZdJw2y_F9jsbD9d7-O_xy/pubhtml?gid=1494122819&single=true)

A large contributing factor to the FP rate is that Manta often has large uncertainty windows around breakpoints in the order of hundreds of bases. This often occurs in regions with multiple structural variants or a high amount of noise (due to repetitive sections), or occasionally for no discernable reason.

BPI uses Manta’s variant calls to re-analyse BAM files and precisely determine the location of the breaks, and applies a set of filters to remove false positives, thereby increasing the accuracy of Manta’s calls.

We also calculate an Allele Frequency which provides an input to determining the ploidy of the variant (when combined with purity).

For convenience, we also output a BAM which has been filtered +/- 500 bp around each breakpoint.
See Also
[Illumina Manta](https://github.com/Illumina/manta)
[Papenfuss Lab Gridss](https://github.com/PapenfussLab/gridss)
[Nature Methods - CONSERTING: integrating copy-number analysis with structural-variation detection](http://www.nature.com/nmeth/journal/v12/n6/full/nmeth.3394.html)

## Method

### Temporary BAM Creation

For every variant, we create a BAM sorted by “QueryName” which allows us to assess both ends of a paired read in a streaming manner.
We will assess secondary/supplementary reads along with their mates as a pair.
Read all Alignments from both NORMAL and TUMOR BAM around the Manta start and end BP confidence interval +/- 500 (default, can be modified via --proximity) bases of each variant
Save to a “QueryName” sorted BAM in a temporary working directory

### Homology Handling

Manta would always take the 5” (relative to reference) position of the first breakpoint, which depending on the orientation of the variant, may or may not include the homology. The convention also varied if the variant was a BND or another type.

For consistency, BPI will always include the homologous sequence in the first breakpoint, and exclude it from the second breakpoint regardless of variant type.

### First Step - Validate and refine breakpoints from TUMOR reads

Internally, BPI will place the breakpoints at both ends closest to where we expect soft-clipping.
We later adjust these locations for homologous sequences according to the above convention.

For a Manta “Precise” Call
1. We will take the breakpoints exactly as Manta has called them.

For a Manta “Imprecise” Call
1. Process all “interesting” paired read evidence
    1 MUST have pair orientation matching expected orientation of the variant
    1 MUST have each read end from the expected chromosome
    1 NOT proper OR has clipping on expected side (relative to variant orientation)
1. Process all SR-only (i.e. normal pairs which are clipped on an expected side)
1. Process all pairs where there is a secondary alignment with clipping on the expected side
1. From all these pairs, we process the clipping and build a list of clipped sequences at each location, and count the reads supporting each
1. From the clipped sequence list, pick the position with the strongest support within the Manta uncertainty window around each break point.
1. In absence of clipping, we pick the inner-most position from “interesting” pairs
1. Report refined breakpoints

### Second Step - Collect Evidence Stats for Tumor & Normal
With accurate breakpoints, we then count:
1. PR (which span the points)
    1. Correct Orientation
    1. The sum of the distance of each alignment to the breakpoint must be <400 bp.
1. SR (which clip the points) evidence
    1. Correct Orientation
    1. Clipping at the precise location of break
1. Filtering is then applied to the counts of (PR Only, PR+SR, SR Only) aggregated at each breakpoint of the SV in addition to other metrics.
1. The Allele Frequency is calculated from these counts

### Filters

A “short variant” is a delete or dupe <1000 bp.

Name | Description
-----|------------
BPI_BreakpointError | BPI failed to determine breakpoints
BPI_MinDepth | The depth across one of the breakpoints is <10
BPI_MinAnchorLength | For short variants, there has to be at least one SR pair which has an alignment with >=30 bp matched. For others, there must be a PR where both alignments >= 30 bp matched.
BPI_SRSupportZero | Short variant must have SR support at BOTH breakpoints
BPI_SRNormalSupport | Short variant has SR support in normal
BPI_PRNormalSupport | PR support in the normal
BPI_PRSupportZero | No PR support in tumor
BPI_ClippingConcordance | At least 5 base clipped bases concordance between tumor and normal -- usually indicative of a mapping artifact

### Allele Frequency

We must determine an AF at each breakpoint, as structural variants occur in addition to other copy number events like whole chromosome duplication or loss.

* For each breakpoint:
    * For Short Deletes / Duplicates:
        * SupportSRNormalSR + SupportSR
    * Otherwise:
        * SupportSR+SupportPRNormalSR+NormalPR+SupportPR + SupportSR

### Known Issues with Manta

* Manta will incorrectly position the break when there is a SNP which disrupts the homology.
* Manta will return the incorrect sequence (seems to be an off-by-one error) for the homology for certain BND orientations and Tandem Duplications.
* Manta varies between 0 basing and 1 basing for it’s breakpoint ends.
* Manta sometimes misses clear SR evidence in the tumor and reports imprecise breakpoints, when the precise breakpoints are clear
* Manta sometimes misses clear SR evidence in the normal which strongly indicate the event is not somatic.

### Future Work

[ ] We will perform a basic local re-assembly in order to validate the variant call.
