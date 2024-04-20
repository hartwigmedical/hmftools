package com.hartwig.hmftools.wisp.purity;

import java.util.Set;

public enum WriteType
{
    CN_DATA,
    CN_PLOT,
    LOH_DATA,
    SOMATIC_DATA,
    SOMATIC_PLOT,
    FRAG_LENGTHS;

    public static final String ALL = "ALL";

    public static boolean plotSomatics(final Set<WriteType> writeTypes) { return writeTypes.contains(SOMATIC_PLOT); }
    public static boolean plotCopyNumber(final Set<WriteType> writeTypes) { return writeTypes.contains(CN_PLOT); }
}
