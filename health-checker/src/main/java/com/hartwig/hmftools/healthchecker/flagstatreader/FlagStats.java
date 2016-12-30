package com.hartwig.hmftools.healthchecker.flagstatreader;

import com.google.common.annotations.VisibleForTesting;

import org.jetbrains.annotations.NotNull;

public class FlagStats {

    @NotNull
    private final FlagStatsType flagStatsType;
    private final double value;

    FlagStats(@NotNull final FlagStatsType flagStatsType, final double value) {
        this.flagStatsType = flagStatsType;
        this.value = value;
    }

    @NotNull
    @VisibleForTesting
    FlagStatsType getFlagStatsType() {
        return flagStatsType;
    }

    public double getValue() {
        return value;
    }
}
