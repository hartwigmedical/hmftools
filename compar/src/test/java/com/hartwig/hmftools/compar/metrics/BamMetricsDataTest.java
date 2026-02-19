package com.hartwig.hmftools.compar.metrics;

import static com.hartwig.hmftools.compar.common.CategoryType.TUMOR_BAM_METRICS;
import static com.hartwig.hmftools.compar.metrics.BamMetricsComparer.coverageString;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.FLD_DUPLICATE_PERCENTAGE;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;

public class BamMetricsDataTest
        extends ComparableItemTest<BamMetricsData, BamMetricsComparer, TestBamMetricsDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new BamMetricsComparer(TUMOR_BAM_METRICS, new ComparConfig());
        builder = TestBamMetricsDataBuilder.BUILDER;
        BamMetricsData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
                FLD_DUPLICATE_PERCENTAGE, b -> b.duplicatePercentage = alternateValueSource.metrics().duplicatePercent(),
                coverageString(30), b -> b.percentage30x = alternateValueSource.metrics().coveragePercent(30),
                coverageString(60), b -> b.percentage60x = alternateValueSource.metrics().coveragePercent(60)
        );
        nameToAlternateIndexInitializer = Collections.emptyMap();
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
    }
}
