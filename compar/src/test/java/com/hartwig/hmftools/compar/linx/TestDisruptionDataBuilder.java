package com.hartwig.hmftools.compar.linx;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.compar.TestComparableItemBuilder;
import com.hartwig.hmftools.compar.common.Category;

public class TestDisruptionDataBuilder
{
    public Category subCategory = Category.DISRUPTION;
    public String geneName = "BRAF";
    public List<BreakendData> breakends = List.of(TestBreakendDataBuilder.BUILDER.create());

    private static final Consumer<TestDisruptionDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.subCategory = Category.DISRUPTION;
        b.geneName = "BRCA2";
        b.breakends = List.of(TestBreakendDataBuilder.BUILDER.create(), TestBreakendDataBuilder.BUILDER.createWithAlternateDefaults());
    };

    public static final TestComparableItemBuilder<TestDisruptionDataBuilder, DisruptionData> BUILDER =
            new TestComparableItemBuilder<>(TestDisruptionDataBuilder::new, TestDisruptionDataBuilder::build, ALTERNATE_INITIALIZER);

    private DisruptionData build()
    {
        return new DisruptionData(subCategory, geneName, breakends);
    }
}
