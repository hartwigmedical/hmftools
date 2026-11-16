package com.hartwig.hmftools.common.metrics;

import java.util.Map;

public class GeneDepth
{
    public final String Gene;
    public final String Chromosome;
    public final int PosStart;
    public final int PosEnd;
    public final Map<DepthRange,Integer> DepthCounts;

    public GeneDepth(
            final String gene, final String chromosome, final int posStart, final int posEnd,
            final Map<DepthRange,Integer> depthCounts)
    {
        Gene = gene;
        Chromosome = chromosome;
        PosStart = posStart;
        PosEnd = posEnd;
        DepthCounts = depthCounts;
    }
}
