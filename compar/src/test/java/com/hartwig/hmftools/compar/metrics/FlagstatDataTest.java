package com.hartwig.hmftools.compar.metrics;

import static com.hartwig.hmftools.compar.common.CategoryType.TUMOR_FLAGSTAT;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.FLD_MAPPED_PROPORTION;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;

public class FlagstatDataTest extends ComparableItemTest<FlagstatData, FlagstatComparer, TestFlagstatDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new FlagstatComparer(TUMOR_FLAGSTAT, new ComparConfig());
        builder = TestFlagstatDataBuilder.BUILDER;
        FlagstatData alternateValueSource = builder.createWithAlternateDefaults();
        fieldToAlternateValueInitializer =
                Map.of(FLD_MAPPED_PROPORTION, b -> b.mappedProportion = alternateValueSource.flagStats().mappedProportion());
        nameToAlternateIndexInitializer = Collections.emptyMap();
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
    }
}
