package com.hartwig.hmftools.virusdetect;

import java.util.Comparator;

// Per-contig support over a set of aligned reads.
public record ContigStats(
        ViralContig contig,
        // Reads with any alignment to this contig
        int readCount,
        // Reads with more than one alignment to this contig (BWA -a repeats/multi-loci)
        int multiAlignReads,
        // Alignments to this contig per read (>= 1)
        SummaryStats alignPerRead,
        // Alignments dropped for clipping over the contig start/end (circular-genome artifact)
        int originClippedReads,
        // Contig positions with at least one aligned base
        int coveredBases,
        // Depth spans the whole contig, so uncovered positions count as depth 0
        SummaryStats depth,
        // BWA alignment score distribution across those reads
        SummaryStats alignerScore,
        // Strain support: reads softly attributed to this strain, split across the sibling contigs by divergence
        double readVotes)
{
    // Best supported first: most read votes, with the contig name breaking ties for determinism.
    public static final Comparator<ContigStats> BEST_SUPPORT_FIRST =
            Comparator.comparingDouble(ContigStats::readVotes).reversed().thenComparing(stats -> stats.contig().name());

    public double coverageFraction()
    {
        int contigLength = contig.length();
        return contigLength == 0 ? 0 : (double) coveredBases / contigLength;
    }
}
