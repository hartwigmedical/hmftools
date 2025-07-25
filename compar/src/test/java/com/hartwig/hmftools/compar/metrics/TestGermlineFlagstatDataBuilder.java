package com.hartwig.hmftools.compar.metrics;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.metrics.ImmutableBamFlagStats;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestGermlineFlagstatDataBuilder
{
    public double mappedProportion = 0.99;

    private static final Consumer<TestGermlineFlagstatDataBuilder> ALTERNATE_INITIALIZER = b -> b.mappedProportion = 0.9;

    public static final TestComparableItemBuilder<TestGermlineFlagstatDataBuilder, GermlineFlagstatData> BUILDER =
            new TestComparableItemBuilder<>(
                    TestGermlineFlagstatDataBuilder::new,
                    TestGermlineFlagstatDataBuilder::build,
                    ALTERNATE_INITIALIZER
            );

    private GermlineFlagstatData build()
    {
        return new GermlineFlagstatData(ImmutableBamFlagStats.builder()
                .mappedProportion(mappedProportion)
                .uniqueReadCount(-1)
                .secondaryCount(-1)
                .supplementaryCount(-1)
                .duplicateProportion(-1)
                .pairedInSequencingProportion(-1)
                .properlyPairedProportion(-1)
                .withItselfAndMateMappedProportion(-1)
                .singletonProportion(-1)
                .build());
    }
}
