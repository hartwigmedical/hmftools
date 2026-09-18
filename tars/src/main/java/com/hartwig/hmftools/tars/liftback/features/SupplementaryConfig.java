package com.hartwig.hmftools.tars.liftback.features;

public class SupplementaryConfig
{
    public final int MinIntronLength;
    public final int MaxIntronLength;
    public final int MaxSuppMerges;
    public final boolean AnnotatedOnly;
    public final int MaxSuppReadOverlap;

    public SupplementaryConfig(
            final int minIntronLength, final int maxIntronLength, final int maxSuppMerges,
            final boolean annotatedOnly, final int maxSuppReadOverlap)
    {
        MinIntronLength = minIntronLength;
        MaxIntronLength = maxIntronLength;
        MaxSuppMerges = maxSuppMerges;
        AnnotatedOnly = annotatedOnly;
        MaxSuppReadOverlap = maxSuppReadOverlap;
    }
}
