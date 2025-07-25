package com.hartwig.hmftools.compar.cider;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.cider.ImmutableCdr3Sequence;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestCiderVdjDataBuilder
{
    public String cdr3Seq = "ACGT";
    public String filter = "PASS";
    public String locus = "IGH";

    private static final Consumer<TestCiderVdjDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.cdr3Seq = "TGTAATGT";
        b.filter = "PARTIAL";
        b.locus = "TRB";
    };

    public static final TestComparableItemBuilder<TestCiderVdjDataBuilder, CiderVdjData> BUILDER =
            new TestComparableItemBuilder<>(TestCiderVdjDataBuilder::new, TestCiderVdjDataBuilder::build, ALTERNATE_INITIALIZER);

    private CiderVdjData build()
    {
        return new CiderVdjData(ImmutableCdr3Sequence.builder()
                .cdr3Seq(cdr3Seq)
                .locus(locus)
                .filter(filter)
                .cdr3AA("")
                .blastnStatus("")
                .minHighQualBaseReads(-1)
                .assignedReads(-1)
                .inFrame(true)
                .containsStop(false)
                .build());
    }
}
