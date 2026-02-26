package com.hartwig.hmftools.compar.isofox;

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
    public String comparisonChromosomeStart = "9";
    public int comparisonPositionStart = 1000;
    public String comparisonChromosomeEnd = "9";
    public int comparisonPositionEnd = 2000;

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
        b.comparisonChromosomeStart = "7";
        b.comparisonPositionStart = 3000;
        b.comparisonChromosomeEnd = "7";
        b.comparisonPositionEnd = 4000;
    };

    public static final TestComparableItemBuilder<TestNovelSpliceJunctionDataBuilder, NovelSpliceJunctionData> BUILDER =
            new TestComparableItemBuilder<>(TestNovelSpliceJunctionDataBuilder::new, TestNovelSpliceJunctionDataBuilder::build, ALTERNATE_INITIALIZER);

    private NovelSpliceJunctionData build()
    {
        final NovelSpliceJunction junction = ImmutableNovelSpliceJunction.builder()
                .geneName(geneName)
                .chromosome(chromosome)
                .junctionStart(junctionStart)
                .junctionEnd(junctionEnd)
                .type(type)
                .fragmentCount(fragmentCount)
                .depthStart(-1)
                .depthEnd(-1)
                .regionStart(regionStart)
                .regionEnd(regionEnd)
                .basesStart("")
                .basesEnd("")
                .cohortFrequency(-1)
                .build();

        return new NovelSpliceJunctionData(junction, new BasePosition(comparisonChromosomeStart, comparisonPositionStart), new BasePosition(comparisonChromosomeEnd, comparisonPositionEnd));
    }
}