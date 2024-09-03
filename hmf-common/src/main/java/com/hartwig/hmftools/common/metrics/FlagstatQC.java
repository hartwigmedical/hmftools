package com.hartwig.hmftools.common.metrics;

import org.jetbrains.annotations.NotNull;

public final class FlagstatQC
{
    public static final double MIN_MAPPED_PROPORTION = 0.95;

    public static boolean pass(@NotNull OldFlagstat flagstat)
    {
        return pass(flagstat.mappedProportion());
    }

    public static boolean pass(double mappedProportion)
    {
        return mappedProportion >= MIN_MAPPED_PROPORTION;
    }
}
