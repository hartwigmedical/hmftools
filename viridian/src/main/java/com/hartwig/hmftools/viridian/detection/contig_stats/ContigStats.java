package com.hartwig.hmftools.viridian.detection.contig_stats;

import com.hartwig.hmftools.viridian.common.SummaryStats;
import com.hartwig.hmftools.viridian.reference.ViralContig;

import org.jetbrains.annotations.Nullable;

// Per-contig alignment statistics.
// Used in two contexts:
//   1. All alignments to all contigs (BWA-MEM -a mode).
//   2. Alignments to just the best contig (after the representative contig selection).
public record ContigStats(
        ViralContig contig,
        // Reads with an alignment to this contig.
        int readCount,
        // Reads with an alignment dropped because it clipped over the contig start/end (circular-genome artifact).
        int originClippedReads,
        // Contig positions aligned by at least 1 alignment.
        int coveredBases,
        // Note uncovered bases count as depth=0.
        SummaryStats depth,
        // BWA alignment score distribution across those reads. Null only when every alignment was origin clipped.
        @Nullable SummaryStats alignerScore)
{
    public ContigStats
    {
        if(readCount < 0)
        {
            throw new IllegalArgumentException("Invalid readCount: " + readCount);
        }
        if(originClippedReads < 0)
        {
            throw new IllegalArgumentException("Invalid originClippedReads: " + originClippedReads);
        }
        if(coveredBases < 0 || coveredBases > contig.length())
        {
            throw new IllegalArgumentException("Invalid coveredBases: " + coveredBases);
        }
    }

    public double coverageFraction()
    {
        return (double) coveredBases / contig.length();
    }
}
