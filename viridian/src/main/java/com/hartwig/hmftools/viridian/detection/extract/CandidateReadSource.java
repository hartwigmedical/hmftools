package com.hartwig.hmftools.viridian.detection.extract;

// Determines which subsets of reads are considered for extraction from the tumor BAM.
public enum CandidateReadSource
{
    // Full unmapped, virus decoy contig, half-mapped, clipped.
    // Identifies most viral reads but is slow because it has to scan all contigs.
    ALL,
    // Fully unmapped, virus decoy contig.
    // Misses reads which may be viral integrations (i.e. split host-virus), but much faster because it avoids scanning all the contigs.
    FULLY_UNMAPPED_AND_VIRUS
}
