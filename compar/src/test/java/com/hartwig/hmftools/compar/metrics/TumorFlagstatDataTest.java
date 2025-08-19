package com.hartwig.hmftools.compar.metrics;

import static com.hartwig.hmftools.compar.metrics.MetricsCommon.FLD_MAPPED_PROPORTION;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;

public class TumorFlagstatDataTest extends ComparableItemTest<TumorFlagstatData, TumorFlagstatComparer, TestTumorFlagstatDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new TumorFlagstatComparer(new ComparConfig());
        builder = TestTumorFlagstatDataBuilder.BUILDER;
        TumorFlagstatData alternateValueSource = builder.createWithAlternateDefaults();
        fieldToAlternateValueInitializer =
                Map.of(FLD_MAPPED_PROPORTION, b -> b.mappedProportion = alternateValueSource.mFlagstat.mappedProportion());
        nameToAlternateIndexInitializer = Collections.emptyMap();
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
    }
}
