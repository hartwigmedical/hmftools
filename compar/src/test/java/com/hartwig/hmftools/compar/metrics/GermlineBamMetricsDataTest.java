package com.hartwig.hmftools.compar.metrics;

import static com.hartwig.hmftools.compar.metrics.GermlineBamMetricsData.FLD_PERCENTAGE_10X;
import static com.hartwig.hmftools.compar.metrics.GermlineBamMetricsData.FLD_PERCENTAGE_20X;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.FLD_DUPLICATE_PERCENTAGE;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;

public class GermlineBamMetricsDataTest
        extends ComparableItemTest<GermlineBamMetricsData, GermlineBamMetricsComparer, TestGermlineBamMetricsDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new GermlineBamMetricsComparer(new ComparConfig());
        builder = TestGermlineBamMetricsDataBuilder.BUILDER;
        GermlineBamMetricsData alternateValueSource = builder.createWithAlternateDefaults();
        fieldToAlternateValueInitializer = Map.of(
                FLD_DUPLICATE_PERCENTAGE, b -> b.duplicatePercentage = alternateValueSource.Metrics.duplicatePercent(),
                FLD_PERCENTAGE_10X, b -> b.percentage10x = alternateValueSource.Metrics.coveragePercent(10),
                FLD_PERCENTAGE_20X, b -> b.percentage20x = alternateValueSource.Metrics.coveragePercent(20)
        );
        nameToAlternateIndexInitializer = Collections.emptyMap();
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
    }
}
