package com.hartwig.hmftools.orange.cohort.percentile;

import com.google.common.collect.Multimap;

import org.jetbrains.annotations.NotNull;

public class PercentileModel {

    @NotNull
    private final Multimap<PercentileType, CohortPercentiles> percentileMap;

    public PercentileModel(@NotNull final Multimap<PercentileType, CohortPercentiles> percentileMap) {
        this.percentileMap = percentileMap;
    }
}
