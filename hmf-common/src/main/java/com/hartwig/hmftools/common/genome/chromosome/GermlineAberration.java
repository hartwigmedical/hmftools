package com.hartwig.hmftools.common.genome.chromosome;

import java.util.Arrays;
import java.util.Set;
import java.util.stream.Collectors;

import org.jetbrains.annotations.NotNull;

public enum GermlineAberration {
    NONE,
    MOSAIC_X,
    KLINEFELTER,
    XYY,
    TRISOMY_X,
    TRISOMY_13,
    TRISOMY_15,
    TRISOMY_18,
    TRISOMY_21;

    @NotNull
    public static String toString(Set<GermlineAberration> aberrations) {
        return aberrations.stream().map(Enum::toString).collect(Collectors.joining(","));
    }

    @NotNull
    public static Set<GermlineAberration> fromString(String line) {
        return Arrays.stream(line.split(",")).map(GermlineAberration::valueOf).collect(Collectors.toSet());
    }

}
