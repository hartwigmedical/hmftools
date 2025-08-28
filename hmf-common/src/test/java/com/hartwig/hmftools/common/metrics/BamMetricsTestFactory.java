package com.hartwig.hmftools.common.metrics;

import java.util.List;

import org.jetbrains.annotations.NotNull;

public final class BamMetricsTestFactory
{
    @NotNull
    public static BamMetricsSummary createMinimalTestWGSMetrics()
    {
        List<Integer> coverageLevels = List.of(10, 20, 30, 60);
        List<Double> coveragePercents = List.of(0.0, 0.0, 0.0, 0.0);

        return ImmutableBamMetricsSummary.builder()
                .totalReads(0)
                .totalRegionBases(0)
                .duplicateReads(0)
                .dualStrandReads(0)
                .meanCoverage(0D)
                .sdCoverage(0D)
                .medianCoverage(0)
                .madCoverage(0)
                .lowMapQualPercent(0D)
                .duplicatePercent(0D)
                .unmappedPercent(0D)
                .lowBaseQualPercent(0D)
                .overlappingReadPercent(0D)
                .cappedCoveragePercent(0D)
                .coverageLevels(coverageLevels)
                .coveragePercents(coveragePercents)
                .build();
    }
}
