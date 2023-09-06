package com.hartwig.hmftools.neo.score;

public class VariantsCounts
{
    public final int Coverage;
    public final int VariantSupport;

    public VariantsCounts(final int coverage, final int variantSupport)
    {
        Coverage = coverage;
        VariantSupport = variantSupport;
    }
}
