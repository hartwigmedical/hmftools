# <ISOFOX>

## Overview
ISOFOX is a tool for counting fragment support for gene and transcript features using genome aligned RNASeq data in tumor samples.   In particular, ISOFOX estimates transcript abundance (including unspliced transcripts) and detects evidence for novel splice junctions and retained introns.    The input for ISOFOX is mapped paired end reads (we use STAR for our aligner).

ISOFOX uses a similar methodology to several previous transcript abundance estimation tools, but may offer several advantages by using a genome based mapping:
* Explicit estimates of the abundance unspliced transcripts in each gene
* Avoids overfitting of 'retained intron' transcripts on the basse
* Individual or combinations of splice junctions which are unique to a transcript will be weighed strongly.  Does not overfit variability of coverage within exons (eg. B2M)

## Install

## Running

### Usage

```
java TODO
```

### Mandatory Arguments

Argument | Description 
---|---
TO DO | 

### Optional Arguments

Argument | Default | Description 
---|---|---
TO DO | |
 
## Algorithm

### 1. Modelling and grouping of spliced and unspliced transcripts

For determining transcript abundance, we consider all transcripts in Ensembl, grouped by gene.    Each gene may have 1 to N transcripts.   Since we use ribosomal depletion to collect RNA we need to  explicitly consider that each gene will have unspliced reads.   Hence, we consider an additional ‘unspliced transcript’ per gene which includes all exonic and intronic segments.   Any fragment that overlaps a region which is intronic on all transcripts is assumed to be unspliced.

Genes that overlap each other on the same chromosome (either sense, anti-sense or shared exons) are considered together as a group so that each fragment which could potentially relate to one of multiple genes is only counted once.     

### 2. Modelling sample specific fragment distribution

The fragment length distibution of the sample is measured by sampling the insert size of up to 1 million genic intronic fragments.   Any fragment with an N in the cigar or which overlaps an exon is excluded from the fragment distribution.  A maximum of 1000 fragments is permitted to be sampled per gene so no individual gene can dominate the sample distribution.   

### 3. Expected rates per 'category' and expected GC distibution per transcript

#### A. Expected rates per catageory

For each transcript in a group of overlapping genes, ISOFOX measures the expected proportion of fragments that have been randomly sampled from that transcript with lengths matching the length distribution of the sample that match a specific subset of transcripts (termed a 'category' in ISOFOX, but generally referred to as an equivalence class in other tools such as Salmon).   For any gene, that contains at least 1 transcript with more than 1 exon an 'UNSPLICED' transcript of that gene is also considered as a indpendent transcript that could be expressed.  

The proportion is calculated by determining which category or set of transcripts that fragments of length {50, 100, 150,200,250,300,350,400,450,500,550} bases starting at each possible base in the transcript in question could be a part of.    This is then weighted by the empirically observed

For example a gene with 2 transcripts (A & B) and an UNSPLICED transcript might have the following expected rates:

Category | 'A' Transcript | 'B' Transcript |'UNSPLICED' Transcript 
---|---|---|---
A Only|0.5|0|0
B only|0|0.2|0
Both A & B|0.1|0.2|0
A & UNSPLICED|0.1|0|0.05
B & UNSPLICED|0|0|0
Both A & B & UNSPLICED|0.3|0.6|0.15
UNSPLICED|0|0|0.8
TOTAL|1.0|1.0|1.0

In this example, 50% of all fragments from transcript A are expected to be uniquely mapped to the A transcript, 30% may be mapped to A,B or UNSPLICED (likely fragments matching a long exon), 10% are expected to be mapped to either A or B but with splicing and the final 10% are expected to be mapped to a region which is either exonic in A or unspliced.     These rates are compared to the observed abundance of each category in subsequent steps to estimate the actual abundance of each transcript.

#### B. Expected GC distribution

For each of the sampled fragments in each transcript (including unspliced transcripts), the GC content is also measured.  For each transcript the unbiased proportion of expected fragments with each GC count between 0 and 2x the read length is calculated.   This is used subsequently to estimate GC bias.

*TODO* - Problem with fragments < 2 x read length.

### 4. Counting abundance per unique group of shared transcripts

Similarly to the estimated rate calculation above we also use the same grouping of transcripts together across all genes which overlap each other to determine actual counts.   We assume that any fragment that overlaps this region must belong either to one of these transcripts or to an unspliced version of one of the genes.

Each fragment is assigned to a 'category' based on the set of transcripts that it may belong to.   We allow a fragment may belong to a transcript if:
* Every base of the fragment is exonic in that transcript (allowing for homology with reads that marginally overlap exon boundaries) AND
* Every splice junction called exists in that transcript AND
* the distance between the paired reads in that transcript is not > maximum insert size distribution 

Any fragment which does not contain a splice junction, is wholly contained within the bounds of a gene, and with fragment size <= maximum insert size distribution is also allowed to map to an ‘UNSPLICED’ transcript of that gene.

Note that reads which marginally overhang an exon boundary or are soft clipped at or beyond an exon boundary have special treatment. This is particularly relevant for reads that have an overhang of 1 or 2 bases which will not be mapped by STAR with default parameters   If the overhanging section can be uniquely mapped either to the reference or to the other side of only a single known spliced junction, then the fragment is deemed to be supporting that splice junction or in the case of supporting just the reference is deemed to be supporting the UNSPLICED transcript.  If it cannot be uniquely mapped or matches neither of those locations exactly, it is truncated at the exon boundary.

We currently ignore marked duplicates in our counting.  This may lead to understimation of abundance for high coverage genes and transcripts which may reach saturation of all read pair positions (particularly for genes where the per base coverage exceeds ~3000). We will address this problem in a later release <TO DO>.

### 5. Fit abundance estimate per transcript

For each group of transcripts conidered together we aim to fit the relative abundance.   Like many previous tools (RSEM, Salmon, Kallisto, etc), we use an expectation maximastion algorithm to find the allocation of fragments to each transcript which give the least residuals compared to the expected rates for each transcript.

<TO DO: Add a step which improves on this by removing or limiting allocation to transcripts where fitted fragments >> observed fragments for private allocations >

### 6. Bias Estimation and Correction

#### A. GC Bias estimate

Expected GC distribution for sample is estimated to be the sum of the distribution for each transcript (calcluated above) multiplied by the proportion of fragments in the sample which have been estimated (nb - this is similar to the methodology used in Salmon).  To avoid overfitting highly expressed genes or transcripts, genes with more than X% of fragments are ignored in this calculation

The GC bias for each integer GC content level between 1 & 2 x read length is the ratio of the estimated to the actual distribution ignoring the highly expressed genes. 

##### B. Fragment Length Bias

TO DO

#### C. Sequence Start Speicfic Bias

TO DO 

#### D. 5' CAP bias

TO DO

#### E. Adjust expected rates for biases and restimate abundances per transcript

TO DO - Apply adjustments to the expected rates.   Then repeat step 5.


### 7. Counting and characterisation of novel splice junctions

A novel splice junction is considered to be a splicing event called by the aligner which is not part of any annotated event.  For each novel splice junction we count the number of fragments supporting the event as well as the total coverage at each end of the splicing junction.    Each novel splice junction is classified as one of the following types of events

* SKIPPED_EXON - both splice sites are splice sites on a single existing transcript but the interim exon(s) are skipped.
* MIXED_TRANSCRIPT - both splice sites are splice sites on different existing transcripts but not the same one
* NOVEL_EXON - One end of the splice junction is cis phased with another novel splice junction to form a novel exon that does not exist on any transcript.
* NOVEL_3_PRIME_SS - 5' end matches a known splice site, but 3' end does not.  3' end may be intronic or exonic
* NOVEL_5_PRIME_SS - 3' end matches a known splice site, but 5' end does not.  5' end may be intronic or exonic
* NOVEL_INTRON -  Neither end matches a known splice site.  Both 5' and 3' ends are wholly contained within a single exon on all transcripts
* INTRONIC -  Neither end matches a known splice site.  Both 5' and 3' ends are intronic
* INTRONIC_TO_EXONIC - Neither end matches a known splice site.  One end is intronic and 1 end is exonic

In the case of overlapping genes, we assign the novel splice junction to one of the genes using the following priority rules in order
* Genes with matching splice site at a least one end
* Genes on strand such that splice motif matches canonical splice motif (GT-AG)
* Gene with most transcripts

For each novel splice junction, we also record the distance to the nearest splice junction for each novel site, the motif of the splice site and the transcripts compatible with either end of the alternative splicing.


### 8. Counting and characterisation of retained introns

We also search explicitly for evidence of retained introns, ie where reads overlap exon boundaries.   We may find many such reads as any unspliced transcripts can have such fragments.    To reduce false positives, we only consider exon boundaries which do not have exons on other transcripts overlapping them, and we hard filter any evidence where we don't see at least 3 reads overlapping the exon boundary or at least 1 read from a fragment that contains another splice junction.

<TO DO - Add filtering that the read count must signiificantly exceed the unspliced coverage of all gene overlapping that base>

### 9. Counting and characterisation of chimeric and read through junctions

<TO DO>

## 10. Annotation of novel splice features

<TO DO - Combine intron retention and novel splice site information>

### 11. Panel of Normals  <TO DO>

We have developed a 'panel of normals' for both novel splice junctions and novel retained introns across a cohort of 1700 samples.  In this context 'novel' means any splice event that does not exist in an ensembl annotated transcript. The panel of normals is created to estimate population level frequencies of each of the 'novel' features.   

For each novel splice junction  we count
* Number of unique samples with novel splice junction
* Total # of fragments supporting novel splice junction across all the unique samples

For intron retention cases we count
* Number of unique samples with at least 3 fragments or at least one fragment with a known splice event supporting intron retention
* Total # of fragments supporting intron retention from those unique samples
* Total # of fragments supporting intron retention from those unique samples which also have splcing.

Each novel splice junction and retained intron for each sample is annotated with the population level frequencies

## Outputs

### Fragment length distribution

Field | Description 
---|---
FragmentLength | 
Count |

### Gene level data

Field | Description 
---|---
Chromosome | 
GeneLength | 
IntronicLength | 
TransCount | 
TotalFragments | 
SupportingTrans | 
Alt | 
Unspliced | 
ReadThrough | 
Chimeric | 
Duplicates | 
UnsplicedAlloc | 
FitResiduals | 
GeneSet | 

### Transcriptome level data


Field | Description 
---|---
ExonCount | 
EffectiveLength | 
ExonsMatched | 
ExonicBases | 
ExonicCoverage | 
FitAllocation | 
UniqueBases | 
UniqueBaseCoverage | 
UniqueBaseAvgDepth | 
SpliceJuncSupported | 
UniqueSpliceJunc | 
UniqueSpliceJuncSupported | 
ShortFragments | 
ShortUniqueFragments | 
LongFragments | 
LongUniqueFragments | 
SpliceJuncFragments | 
UniqueSpliceJuncFragments | 

### Expected abundance rates per transcript per category

Field | Description 
---|---
Transcript |
Category |
Rate |

### Observed and fitted counts per category

Field | Description 
---|---
Category |
Count |
FittedCount |

### Novel splice junctions

Field | Description 
---|---
Chromosome | 
Strand | 
SjStart | 
SjEnd | 
FragCount | 
StartDepth | 
EndDepth | 
Type | 
StartContext | 
EndContext | 
NearestStartExon | 
NearestEndExon | 
StartBases | 
EndBases | 
StartTrans | 
EndTrans | 


### Novel retained introns

Field | Description 
---|---
Chromosome|
Strand|
Position|
Type|
FragCount|
SplicedFragCount|
TotalDepth|
TranscriptInfo|
