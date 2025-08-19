package com.hartwig.hmftools.compar.cider;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.cider.ImmutableCdr3LocusSummary;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestCdr3LocusSummaryDataBuilder
{
    public String locus = "IGH";
    public int passSequenceCount = 10000;

    private static final Consumer<TestCdr3LocusSummaryDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.locus = "TRB";
        b.passSequenceCount = 13000;
    };

    public static final TestComparableItemBuilder<TestCdr3LocusSummaryDataBuilder, Cdr3LocusSummaryData> BUILDER =
            new TestComparableItemBuilder<>(
                    TestCdr3LocusSummaryDataBuilder::new,
                    TestCdr3LocusSummaryDataBuilder::build,
                    ALTERNATE_INITIALIZER
            );

    private Cdr3LocusSummaryData build()
    {
        return new Cdr3LocusSummaryData(ImmutableCdr3LocusSummary.builder()
                .locus(locus)
                .passSequences(passSequenceCount)
                .readsUsed(-1)
                .readsTotal(-1)
                .downSampled(false)
                .sequences(-1)
                .build());
    }
}
