package com.hartwig.hmftools.virusdetect;

// Per-contig support over a set of aligned reads, plus how it fares against the sibling contigs of the same virus.
public record ContigStats(
        String contig,
        int contigLength,
        int readCount,              // reads with any alignment to this contig
        int coveredBases,           // contig positions with at least one aligned base
        // Depth spans the whole contig, so uncovered positions count as depth 0.
        int minDepth,
        int maxDepth,
        double meanDepth,
        double meanAlignerScore,    // mean BWA alignment score across those reads
        // Strain rivalry against the near-identical sibling contigs of the same virus:
        double readVotes,           // reads softly attributed to this strain, split across contigs by divergence
        int readsBestInRivals,      // reads won outright where a rival contig was also in play
        // Divergence lead over the runner-up on the won reads; NaN when no read was contested.
        double marginMean,
        double marginMedian,
        double marginP90)
{
    public double coverageFraction()
    {
        return contigLength == 0 ? 0 : (double) coveredBases / contigLength;
    }
}
