package com.hartwig.hmftools.virusdetect;

// Per-contig support over a set of aligned reads. Depth statistics span the entire contig, so uncovered
// positions count as depth 0.
public record ContigStats(
        String contig, int contigLength, int readCount, int coveredBases,
        int minDepth, int maxDepth, double meanDepth, double meanAlignerScore)
{
    public double coverageFraction()
    {
        return contigLength == 0 ? 0 : (double) coveredBases / contigLength;
    }
}
