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

### Mandatory Arguments

### Optional Arguments

## Algorithm

### Modelling spliced and unspliced transcripts

For determining transcript abundance, we consider all transcripts in Ensembl, grouped by gene.    Each gene may have 1 to N transcripts.   Since we use ribosomal depletion to collect RNA we need to  explicitly consider that each gene will have unspliced reads.   Hence, we consider an additional ‘unspliced transcript’ per gene which includes all exonic and intronic segments.   Any fragment that overlaps a region which is intronic on all transcripts is assumed to be unspliced.

Genes that overlap each other(either sense, anti-sense or shared exons) are considered together as a group so that each fragment is only counted once.     

<TO DO: NOTE on handling of duplicates>

### Modelling sample specific fragment distribution

<Using exclusively reads that overlap exon.    Exclude reads with N in cigar and that unusual gene with huge numbers>



### Calculate expected shared and private abundance rates per transcript


Counting abundance per unique group of shared transcripts
For counting fragments, we first group transcripts together across all genes which overlap each other at all in either in intronic or exonic regions.   Any fragment that overlaps this region must belong either to one of these transcripts or to an unspliced version of one of the genes.

We consider that a fragment may belong to a transcript if:
Every base of the fragment is exonic in that transcript (allowing for homology with reads that marginally overlap exon boundaries) AND
Every splice junction called exists in that transcript AND
the distance between the paired reads in that transcript is not > maximum insert size distribution 

Any fragment which does not contain a splice junction, is wholly contained within the bounds of a gene, and with fragment size <= maximum insert size distribution is also allowed to map to an ‘UNSPLICED’ transcript of that gene.

Each fragment is assigned to a category 

### Fit abundance estimate per transcript

Like many previous tools (RSEM, Salmon, Kallisto, etc), we have chosen to use an expectations maximisation algorithm . . .

### Bias Estimation and Correction

<TO DO: GC Bias, Fragment Length Bias, Sequence Start Speicfic Bias, 5' CAP bias>

### Counting and characterisation of novel splice junctions

TO DO

### Counting and characterisation of chimeric and read through junctions

## Annotation of novel splice features

### Panel of Normals  (TO DO)

We have developed a 'panel of normals' for both novel splice junctions and novel retained introns across a cohort of 1700 samples.  In this context 'novel' means any splice event that does not exist in an ensembl annotated transcript. The panel of normals is created to estimate population level frequencies of each of the 'novel' features.   

For each novel splice junction  we count
* # of unique samples with novel splice junction
* Total # of fragments supporting novel splice junction across all the unique samples

For intron retention cases we count
* # of unique samples with at least 3 fragments or at least one fragment with a known splice event supporting intron retention
* Total # of fragments supporting intron retention from those unique samples
* Total # of fragments supporting intron retention from those unique samples which also have splcing.

Each novel splice junction and retained intron for each sample is annotated with the population level frequencies

## Outputs

### Fragment length distribution

### Gene level data

### Transcriptome level data

### Expected abundance rates per transcript per category

### Observed and fitted counts per category

### Novel splice junctions

### Novel retained introns

