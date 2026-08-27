package com.hartwig.hmftools.compar.cuppa;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;

public class CuppaDataTest extends ComparableItemTest<CuppaData, CuppaComparer, TestCuppaDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new CuppaComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestCuppaDataBuilder.BUILDER;
        CuppaData alternateValueSource = builder.createWithAlternateDefaults();
        fieldToAlternateValueInitializer = Map.of(
                CuppaComparer.Fields.TopCancerType.toString(), b -> b.cancerType = alternateValueSource.PredictionEntry.CancerType,
                CuppaComparer.Fields.Probability.toString(), b -> b.dataValue = alternateValueSource.PredictionEntry.DataValue
        );
        nameToAlternateIndexInitializer = Map.of(
                "dataType", b -> b.dataType = alternateValueSource.PredictionEntry.DataType,
                "classifierName", b -> b.classifierName = alternateValueSource.PredictionEntry.ClassifierName
        );
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
    }
}
