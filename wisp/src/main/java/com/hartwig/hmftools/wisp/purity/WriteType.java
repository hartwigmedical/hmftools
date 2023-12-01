package com.hartwig.hmftools.wisp.purity;

import java.util.Set;

public enum WriteType
{
    CN_DATA,
    CN_PLOTS,
    SOMATICS,
    SOMATIC_ALL,
    SOMATIC_PEAK,
    SOMATIC_PLOT;

    public static final String ALL = "ALL";

    public static boolean plotSomatics(final Set<WriteType> writeTypes)
    {
        return writeTypes.contains(SOMATIC_PLOT) ||  writeTypes.contains(SOMATIC_ALL);
    }

    public static boolean plotCopyNumber(final Set<WriteType> writeTypes)
    {
        return writeTypes.contains(CN_PLOTS) && writeTypes.contains(CN_DATA);
    }
}
