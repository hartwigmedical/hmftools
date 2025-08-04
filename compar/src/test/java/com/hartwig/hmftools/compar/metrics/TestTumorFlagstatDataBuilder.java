package com.hartwig.hmftools.compar.metrics;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.metrics.ImmutableBamFlagStats;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestTumorFlagstatDataBuilder
{
    public double mappedProportion = 0.99;

    private static final Consumer<TestTumorFlagstatDataBuilder> ALTERNATE_INITIALIZER = b -> b.mappedProportion = 0.9;

    public static final TestComparableItemBuilder<TestTumorFlagstatDataBuilder, TumorFlagstatData> BUILDER =
            new TestComparableItemBuilder<>(TestTumorFlagstatDataBuilder::new, TestTumorFlagstatDataBuilder::build, ALTERNATE_INITIALIZER);

    private TumorFlagstatData build()
    {
        return new TumorFlagstatData(ImmutableBamFlagStats.builder()
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
