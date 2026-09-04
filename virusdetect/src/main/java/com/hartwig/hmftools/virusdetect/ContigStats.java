package com.hartwig.hmftools.virusdetect;

// Per-contig support over a set of aligned reads. Depth spans the whole contig, so uncovered positions count as
// depth 0. The vote and margin fields quantify strain rivalry between the near-identical contigs of one virus: how
// many reads probabilistically belong to this strain, and, over reads contested between strains, how decisively it
// out-aligns the runner-up. Margins are absent when no read contests the contig.
public record ContigStats(
        String contig, int contigLength, int readCount, int coveredBases,
        int minDepth, int maxDepth, double meanDepth, double meanAlignerScore,
        double voteReads, int contestedReads, double marginMean, double marginMedian, double marginP90)
{
    public double coverageFraction()
    {
        return contigLength == 0 ? 0 : (double) coveredBases / contigLength;
    }
}
