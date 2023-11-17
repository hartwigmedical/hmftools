package com.hartwig.hmftools.wisp.purity;

import java.util.Set;

public enum WriteType
{
    CN_DATA,
    CN_PLOTS,
    SOMATICS,
    SOMATICS_ALL,
    SOMATIC_PLOTS;

    public static final String ALL = "ALL";

    public static boolean plotSomatics(final Set<WriteType> writeTypes)
    {
        return writeTypes.contains(SOMATIC_PLOTS) && (writeTypes.contains(SOMATICS) || writeTypes.contains(SOMATICS_ALL));
    }

    public static boolean plotCopyNumber(final Set<WriteType> writeTypes)
    {
        return writeTypes.contains(CN_PLOTS) && writeTypes.contains(CN_DATA);
    }
}
