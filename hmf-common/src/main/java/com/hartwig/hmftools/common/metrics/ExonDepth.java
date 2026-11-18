package com.hartwig.hmftools.common.metrics;

import java.util.Map;

public class ExonDepth extends GeneDepth
{
    public final int ExonRank;
    public final double MedianDepth;
    public final double MeanDepth;

    public ExonDepth(
            final String gene, final String chromosome, final int posStart, final int posEnd, final int exonRank,
            final double medianDepth, final double meanDepth, final Map<DepthRange,Integer> depthCounts)
    {
        super(gene, chromosome, posStart, posEnd, depthCounts);
        ExonRank = exonRank;
        MedianDepth = medianDepth;
        MeanDepth = meanDepth;
    }
}
