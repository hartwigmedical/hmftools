package com.hartwig.hmftools.compar.metrics;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.metrics.BamFlagStats;
import com.hartwig.hmftools.common.metrics.ImmutableBamFlagStats;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;
import com.hartwig.hmftools.compar.common.CategoryType;

public class TestFlagstatDataBuilder
{
    public double mappedProportion = 0.99;

    private static final Consumer<TestFlagstatDataBuilder> ALTERNATE_INITIALIZER = b -> b.mappedProportion = 0.9;

    public static final TestComparableItemBuilder<TestFlagstatDataBuilder, FlagstatData> BUILDER =
            new TestComparableItemBuilder<>(TestFlagstatDataBuilder::new, TestFlagstatDataBuilder::build, ALTERNATE_INITIALIZER);

    private FlagstatData build()
    {
        BamFlagStats bamFlagStats = ImmutableBamFlagStats.builder()
                .mappedProportion(mappedProportion)
                .uniqueReadCount(-1)
                .secondaryCount(-1)
                .supplementaryCount(-1)
                .duplicateProportion(-1)
                .pairedInSequencingProportion(-1)
                .properlyPairedProportion(-1)
                .withItselfAndMateMappedProportion(-1)
                .singletonProportion(-1)
                .build();

        return new FlagstatData(CategoryType.TUMOR_FLAGSTAT, bamFlagStats);
    }
}
