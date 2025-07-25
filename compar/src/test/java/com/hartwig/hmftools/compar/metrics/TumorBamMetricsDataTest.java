package com.hartwig.hmftools.compar.metrics;

import static com.hartwig.hmftools.compar.metrics.MetricsCommon.FLD_DUPLICATE_PERCENTAGE;
import static com.hartwig.hmftools.compar.metrics.TumorBamMetricsData.FLD_PERCENTAGE_30X;
import static com.hartwig.hmftools.compar.metrics.TumorBamMetricsData.FLD_PERCENTAGE_60X;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;

public class TumorBamMetricsDataTest
        extends ComparableItemTest<TumorBamMetricsData, TumorBamMetricsComparer, TestTumorBamMetricsDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new TumorBamMetricsComparer(new ComparConfig());
        builder = TestTumorBamMetricsDataBuilder.BUILDER;
        TumorBamMetricsData alternateValueSource = builder.createWithAlternateDefaults();
        fieldToAlternateValueInitializer = Map.of(
                FLD_DUPLICATE_PERCENTAGE, b -> b.duplicatePercentage = alternateValueSource.Metrics.duplicatePercent(),
                FLD_PERCENTAGE_30X, b -> b.percentage30x = alternateValueSource.Metrics.coveragePercent(30),
                FLD_PERCENTAGE_60X, b -> b.percentage60x = alternateValueSource.Metrics.coveragePercent(60)
        );
        nameToAlternateIndexInitializer = Collections.emptyMap();
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
    }
}
