package com.hartwig.hmftools.compar.metrics;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.metrics.ImmutableBamMetricsSummary;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestGermlineBamMetricsDataBuilder
{
    public double duplicatePercentage = 0.1;
    public double percentage10x = 0.95;
    public double percentage20x = 0.9;

    private static final Consumer<TestGermlineBamMetricsDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.duplicatePercentage = 0.3;
        b.percentage10x = 0.8;
        b.percentage20x = 0.75;
    };

    public static final TestComparableItemBuilder<TestGermlineBamMetricsDataBuilder, GermlineBamMetricsData> BUILDER =
            new TestComparableItemBuilder<>(
                    TestGermlineBamMetricsDataBuilder::new,
                    TestGermlineBamMetricsDataBuilder::build,
                    ALTERNATE_INITIALIZER
            );

    private GermlineBamMetricsData build()
    {
        return new GermlineBamMetricsData(ImmutableBamMetricsSummary.builder()
                .duplicatePercent(duplicatePercentage)
                .coverageLevels(List.of(10, 20))
                .coveragePercents(List.of(percentage10x, percentage20x))
                .totalRegionBases(-1)
                .totalReads(-1)
                .duplicateReads(-1)
                .dualStrandReads(-1)
                .meanCoverage(-1)
                .sdCoverage(-1)
                .medianCoverage(-1)
                .madCoverage(-1)
                .lowMapQualPercent(-1)
                .unmappedPercent(-1)
                .lowBaseQualPercent(-1)
                .overlappingReadPercent(-1)
                .cappedCoveragePercent(-1)
                .build());
    }
}
