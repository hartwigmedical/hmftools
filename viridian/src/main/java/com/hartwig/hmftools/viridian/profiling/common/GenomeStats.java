package com.hartwig.hmftools.viridian.profiling.common;

import com.hartwig.hmftools.viridian.reference.VirusGenome;

import org.jetbrains.annotations.Nullable;

// Per-virus-genome alignment statistics.
// Used in two contexts:
//   1. All alignments to all virus genomes (BWA-MEM -a mode).
//   2. Alignments to just the best genome (after representative selection).
public record GenomeStats(
        VirusGenome genome,
        // Reads with an alignment to this virus genome.
        int readCount,
        // Reads with an alignment dropped because they clipped over the virus genome start/end (circular-genome artifact).
        int originClippedReads,
        // Virus genome positions aligned by at least 1 alignment.
        int coveredBases,
        // Note uncovered bases count as depth=0.
        SummaryStats depth,
        // BWA alignment score distribution across those reads. Null only when every alignment was origin clipped.
        @Nullable SummaryStats alignerScore)
{
    public GenomeStats
    {
        if(readCount < 0)
        {
            throw new IllegalArgumentException("Invalid readCount: " + readCount);
        }
        if(originClippedReads < 0)
        {
            throw new IllegalArgumentException("Invalid originClippedReads: " + originClippedReads);
        }
        if(coveredBases < 0 || coveredBases > genome.length())
        {
            throw new IllegalArgumentException("Invalid coveredBases: " + coveredBases);
        }
    }

    public double coverageFraction()
    {
        return (double) coveredBases / genome.length();
    }
}
