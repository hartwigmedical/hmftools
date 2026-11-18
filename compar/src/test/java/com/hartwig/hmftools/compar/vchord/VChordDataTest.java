package com.hartwig.hmftools.compar.vchord;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;

public class VChordDataTest extends ComparableItemTest<VChordData, VChordComparer, TestVChordDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new VChordComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestVChordDataBuilder.BUILDER;
        VChordData alternateValueSource = builder.createWithAlternateDefaults();
        fieldToAlternateValueInitializer = Map.of(
                VChordComparer.Fields.BreastCancerHrdScore.toString(), b -> b.breast = alternateValueSource.Prediction.breastCancerHrdScore(),
                VChordComparer.Fields.OvarianCancerHrdScore.toString(), b -> b.ovarian = alternateValueSource.Prediction.ovarianCancerHrdScore(),
                VChordComparer.Fields.OtherCancerScore.toString(), b -> b.other = alternateValueSource.Prediction.otherCancerScore()
        );
        nameToAlternateIndexInitializer = Collections.emptyMap();
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
    }
}
