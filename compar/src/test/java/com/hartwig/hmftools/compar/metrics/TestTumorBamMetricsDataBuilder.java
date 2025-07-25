package com.hartwig.hmftools.compar.metrics;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.metrics.ImmutableBamMetricsSummary;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestTumorBamMetricsDataBuilder
{
    public double duplicatePercentage = 0.1;
    public double percentage30x = 0.95;
    public double percentage60x = 0.9;

    private static final Consumer<TestTumorBamMetricsDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.duplicatePercentage = 0.3;
        b.percentage30x = 0.8;
        b.percentage60x = 0.75;
    };

    public static final TestComparableItemBuilder<TestTumorBamMetricsDataBuilder, TumorBamMetricsData> BUILDER =
            new TestComparableItemBuilder<>(
                    TestTumorBamMetricsDataBuilder::new,
                    TestTumorBamMetricsDataBuilder::build,
                    ALTERNATE_INITIALIZER
            );

    private TumorBamMetricsData build()
    {
        return new TumorBamMetricsData(ImmutableBamMetricsSummary.builder()
                .duplicatePercent(duplicatePercentage)
                .coverageLevels(List.of(30, 60))
                .coveragePercents(List.of(percentage30x, percentage60x))
                .totalRegionBases(-1)
                .totalReads(-1)
                .duplicateReads(-1)
                .dualStrandReads(-1)
                .meanCoverage(-1)
                .sdCoverage(-1)
                .medianCoverage(-1)
                .madCoverage(-1)
                .lowMapQualPercent(-1)
                .unpairedPercent(-1)
                .lowBaseQualPercent(-1)
                .overlappingReadPercent(-1)
                .cappedCoveragePercent(-1)
                .build());
    }
}
