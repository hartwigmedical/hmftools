package com.hartwig.hmftools.compar.metrics;

import static com.hartwig.hmftools.compar.common.CategoryType.TUMOR_BAM_METRICS;
import static com.hartwig.hmftools.compar.metrics.BamMetricsComparer.TUMOR_COVERAGE_PERCENTAGES;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.metrics.BamMetricSummary;
import com.hartwig.hmftools.common.metrics.ImmutableBamMetricSummary;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestBamMetricsDataBuilder
{
    public double duplicatePercentage = 0.1;
    public double percentage30x = 0.95;
    public double percentage60x = 0.9;

    private static final Consumer<TestBamMetricsDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.duplicatePercentage = 0.3;
        b.percentage30x = 0.8;
        b.percentage60x = 0.75;
    };

    public static final TestComparableItemBuilder<TestBamMetricsDataBuilder, BamMetricsData> BUILDER =
            new TestComparableItemBuilder<>(
                    TestBamMetricsDataBuilder::new,
                    TestBamMetricsDataBuilder::build,
                    ALTERNATE_INITIALIZER
            );

    private BamMetricsData build()
    {
        BamMetricSummary bamMetricSummary = ImmutableBamMetricSummary.builder()
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
                .unmappedPercent(-1)
                .lowBaseQualPercent(-1)
                .overlappingReadPercent(-1)
                .cappedCoveragePercent(-1)
                .build();

        return new BamMetricsData(TUMOR_BAM_METRICS, bamMetricSummary, TUMOR_COVERAGE_PERCENTAGES);
    }
}
