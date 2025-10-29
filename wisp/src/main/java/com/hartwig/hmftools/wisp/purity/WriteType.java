package com.hartwig.hmftools.wisp.purity;

import java.util.List;
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

    public static final List<WriteType> SOMATIC_WRITE_TYPES = List.of(SOMATIC_DATA, SOMATIC_PLOT, FRAG_LENGTHS);
    public static final List<WriteType> LOH_WRITE_TYPES = List.of(LOH_DATA);
    public static final List<WriteType> CN_WRITE_TYPES = List.of(CN_DATA, CN_PLOT);

    public static boolean plotSomatics(final Set<WriteType> writeTypes) { return writeTypes.contains(SOMATIC_PLOT); }
    public static boolean plotCopyNumber(final Set<WriteType> writeTypes) { return writeTypes.contains(CN_PLOT); }
}
