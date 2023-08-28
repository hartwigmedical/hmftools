package com.hartwig.hmftools.common.genome.chromosome;

import java.util.Arrays;
import java.util.Set;
import java.util.stream.Collectors;

public enum GermlineAberration
{
    NONE,
    MOSAIC_X,
    KLINEFELTER,
    XYY,
    TRISOMY_X,
    TRISOMY_13,
    TRISOMY_15,
    TRISOMY_18,
    TRISOMY_21,
    TERTRASOMY_9;

    public static String toString(final Set<GermlineAberration> aberrations)
    {
        return aberrations.stream().map(Enum::toString).collect(Collectors.joining(","));
    }

    public static Set<GermlineAberration> fromString(final String line)
    {
        return Arrays.stream(line.split(",")).map(GermlineAberration::valueOf).collect(Collectors.toSet());
    }
}
