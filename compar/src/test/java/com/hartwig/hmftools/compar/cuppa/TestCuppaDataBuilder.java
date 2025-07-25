package com.hartwig.hmftools.compar.cuppa;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.cuppa.ClassifierGroup;
import com.hartwig.hmftools.common.cuppa.ClassifierName;
import com.hartwig.hmftools.common.cuppa.CuppaPredictionEntry;
import com.hartwig.hmftools.common.cuppa.DataType;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestCuppaDataBuilder
{
    public DataType dataType = DataType.PROB;
    public ClassifierName classifierName = ClassifierName.DNA_COMBINED;
    public String cancerType = "HPB: Liver";
    public double dataValue = 0.9;

    private static final Consumer<TestCuppaDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.dataType = DataType.FEAT_CONTRIB;
        b.classifierName = ClassifierName.COMBINED;
        b.cancerType = "Colorectum/Small intestine/Appendix";
        b.dataValue = 0.1;
    };

    public static final TestComparableItemBuilder<TestCuppaDataBuilder, CuppaData> BUILDER =
            new TestComparableItemBuilder<>(TestCuppaDataBuilder::new, TestCuppaDataBuilder::build, ALTERNATE_INITIALIZER);

    private CuppaData build()
    {
        return new CuppaData(new CuppaPredictionEntry(
                "",
                dataType,
                ClassifierGroup.COMBINED,
                classifierName,
                "",
                -1.,
                cancerType,
                dataValue,
                0,
                0
        ));
    }
}
