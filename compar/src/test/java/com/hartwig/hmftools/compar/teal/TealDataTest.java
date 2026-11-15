package com.hartwig.hmftools.compar.teal;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;

public class TealDataTest extends ComparableItemTest<TealData, TealComparer, TestTealDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new TealComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestTealDataBuilder.BUILDER;
        TealData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer =
                Map.of(TealComparer.Fields.TelomereLength.toString(), b -> b.telomereLength = alternateValueSource.TelomereLength.finalTelomereLength());
        nameToAlternateIndexInitializer = Map.of("Type", b -> b.type = alternateValueSource.TelomereLength.type());
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
    }
}
