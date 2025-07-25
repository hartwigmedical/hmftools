package com.hartwig.hmftools.compar.teal;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.teal.ImmutableTelomereLength;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestTealDataBuilder
{
    public String type = "tumor";
    public double telomereLength = 10000;

    private static final Consumer<TestTealDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.type = "ref";
        b.telomereLength = 30000;
    };

    public static final TestComparableItemBuilder<TestTealDataBuilder, TealData> BUILDER =
            new TestComparableItemBuilder<>(TestTealDataBuilder::new, TestTealDataBuilder::build, ALTERNATE_INITIALIZER);

    private TealData build()
    {
        return new TealData(
                ImmutableTelomereLength.builder()
                        .type(type)
                        .finalTelomereLength(telomereLength)
                        .rawTelomereLength(-1)
                        .fullFragments(-1)
                        .cRichPartialFragments(-1)
                        .gRichPartialFragments(-1)
                        .totalTelomericReads(-1)
                        .purity(-1)
                        .ploidy(-1)
                        .duplicateProportion(-1)
                        .meanReadDepth(-1)
                        .gc50ReadDepth(-1)
                        .build()
        );
    }
}
