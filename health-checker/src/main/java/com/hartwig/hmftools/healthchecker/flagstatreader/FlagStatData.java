package com.hartwig.hmftools.healthchecker.flagstatreader;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;

import org.jetbrains.annotations.NotNull;

public class FlagStatData {

    @NotNull
    private final List<FlagStats> passedStats;
    @NotNull
    private final List<FlagStats> failedStats;

    FlagStatData(@NotNull final List<FlagStats> passedStats, @NotNull final List<FlagStats> failedStats) {
        this.passedStats = passedStats;
        this.failedStats = failedStats;
    }

    @NotNull
    public List<FlagStats> getPassedStats() {
        return passedStats;
    }

    @NotNull
    @VisibleForTesting
    List<FlagStats> getFailedStats() {
        return failedStats;
    }
}
