package com.hartwig.hmftools.compar.isofox;

import java.util.Collections;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionContext;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionType;
import com.hartwig.hmftools.common.rna.ImmutableNovelSpliceJunction;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestNovelSpliceJunctionDataBuilder
{
    public String geneName = "CDKN2A";
    public String chromosome = "9";
    public int junctionStart = 1000;
    public int junctionEnd = 2000;
    public AltSpliceJunctionType type = AltSpliceJunctionType.SKIPPED_EXONS;
    public int fragmentCount = 1000;
    public AltSpliceJunctionContext regionStart = AltSpliceJunctionContext.SPLICE_JUNC;
    public AltSpliceJunctionContext regionEnd = AltSpliceJunctionContext.SPLICE_JUNC;

    private static final Consumer<TestNovelSpliceJunctionDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.geneName = "BRAF";
        b.chromosome = "7";
        b.junctionStart = 3000;
        b.junctionEnd = 4000;
        b.type = AltSpliceJunctionType.NOVEL_5_PRIME;
        b.fragmentCount = 500;
        b.regionStart = AltSpliceJunctionContext.EXONIC;
        b.regionEnd = AltSpliceJunctionContext.INTRONIC;
    };

    public static final TestComparableItemBuilder<TestNovelSpliceJunctionDataBuilder, RnaNovelSpliceJunctionData> BUILDER =
            new TestComparableItemBuilder<>(TestNovelSpliceJunctionDataBuilder::new, TestNovelSpliceJunctionDataBuilder::build, ALTERNATE_INITIALIZER);

    private RnaNovelSpliceJunctionData build()
    {
        final NovelSpliceJunction junction = ImmutableNovelSpliceJunction.builder()
                .geneName(geneName)
                .chromosome(chromosome)
                .junctionStart(junctionStart)
                .junctionEnd(junctionEnd)
                .type(type)
                .transcriptStart("")
                .transcriptEnd("")
                .exonStart(1)
                .exonEnd(2)
                .fragmentCount(fragmentCount)
                .depthStart(-1)
                .depthEnd(-1)
                .regionStart(regionStart)
                .regionEnd(regionEnd)
                .basesStart("")
                .basesEnd("")
                .cohortFrequency(-1)
                .build();

        return RnaNovelSpliceJunctionData.from(
                junction,
                new RnaNovelSpliceJunctionComparer(null, Collections.emptyMap()).fieldsList());
    }
}
