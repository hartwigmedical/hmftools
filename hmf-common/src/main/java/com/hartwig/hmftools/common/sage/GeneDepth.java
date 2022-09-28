package com.hartwig.hmftools.common.sage;

public class GeneDepth
{
    public final String Gene;
    public final String Chromosome;
    public final int PosStart;
    public final int PosEnd;
    public final double MissedVariantLikelihood;
    public final int[] DepthCounts;

    public GeneDepth(
            final String gene, final String chromosome, final int posStart, final int posEnd,
            final double missedVariantLikelihood, final int[] depthCounts)
    {
        Gene = gene;
        Chromosome = chromosome;
        PosStart = posStart;
        PosEnd = posEnd;
        MissedVariantLikelihood = missedVariantLikelihood;
        DepthCounts = depthCounts;
    }
}
