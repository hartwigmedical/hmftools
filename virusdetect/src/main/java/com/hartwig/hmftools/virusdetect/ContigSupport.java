package com.hartwig.hmftools.virusdetect;

import org.jetbrains.annotations.Nullable;

// Per-contig support information from alignment.
public record ContigSupport(
        ViralContig contig,
        ContigFilterStatus filterStatus,
        // Reads with any alignment to this contig.
        int readCount,
        // Reads with more than one alignment to this contig (BWA -a repeats/multi-loci).
        int multiAlignReads,
        // Alignments to this contig per read (>= 1). Null when no read was retained here.
        @Nullable SummaryStats alignPerRead,
        // Reads with an alignment dropped for clipping over the contig start/end (circular-genome artifact).
        int originClippedReads,
        // Contig positions aligned by at least 1 alignment.
        int coveredBases,
        // Note uncovered bases count as depth=0.
        SummaryStats depth,
        // BWA alignment score distribution across those reads. Null when no read was retained here.
        @Nullable SummaryStats alignerScore,
        // Read attribution taking into account alignment edit distance (divergence).
        double readVotes)
{
    public boolean isCandidate()
    {
        return filterStatus == ContigFilterStatus.CANDIDATE;
    }

    public double coverageFraction()
    {
        return coverageFraction(coveredBases, contig);
    }

    public static double coverageFraction(int coveredBases, ViralContig contig)
    {
        return contig.length() == 0 ? 0 : (double) coveredBases / contig.length();
    }
}
