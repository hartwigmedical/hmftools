package com.hartwig.hmftools.compar.teal;

import static com.hartwig.hmftools.compar.teal.TealData.FLD_TELOMERE_LENGTH;
import static com.hartwig.hmftools.compar.teal.TealData.FLD_TYPE;

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
        comparer = new TealComparer(new ComparConfig());
        builder = TestTealDataBuilder.BUILDER;
        TealData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer =
                Map.of(FLD_TELOMERE_LENGTH, b -> b.telomereLength = alternateValueSource.TelomereLength.finalTelomereLength());
        nameToAlternateIndexInitializer = Map.of(FLD_TYPE, b -> b.type = alternateValueSource.TelomereLength.type());
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
    }
}
