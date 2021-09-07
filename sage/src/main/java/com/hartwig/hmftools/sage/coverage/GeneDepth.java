package com.hartwig.hmftools.sage.coverage;

public class GeneDepth
{
    public final String Gene;
    public final double MissedVariantLikelihood;
    public final int[] DepthCounts;

    public GeneDepth(final String gene, final double missedVariantLikelihood, final int[] depthCounts)
    {
        Gene = gene;
        MissedVariantLikelihood = missedVariantLikelihood;
        DepthCounts = depthCounts;
    }

    public String gene() { return Gene; }
}
