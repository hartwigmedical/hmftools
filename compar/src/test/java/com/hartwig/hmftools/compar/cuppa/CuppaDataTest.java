package com.hartwig.hmftools.compar.cuppa;

import static com.hartwig.hmftools.compar.cuppa.CuppaData.FLD_PROBABILITY;
import static com.hartwig.hmftools.compar.cuppa.CuppaData.FLD_TOP_CANCER_TYPE;

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
        comparer = new CuppaComparer(new ComparConfig());
        builder = TestCuppaDataBuilder.BUILDER;
        CuppaData alternateValueSource = builder.createWithAlternateDefaults();
        fieldToAlternateValueInitializer = Map.of(
                FLD_TOP_CANCER_TYPE, b -> b.cancerType = alternateValueSource.PredictionEntry.CancerType,
                FLD_PROBABILITY, b -> b.dataValue = alternateValueSource.PredictionEntry.DataValue
        );
        nameToAlternateIndexInitializer = Map.of(
                "dataType", b -> b.dataType = alternateValueSource.PredictionEntry.DataType,
                "classifierName", b -> b.classifierName = alternateValueSource.PredictionEntry.ClassifierName
        );
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
    }
}
