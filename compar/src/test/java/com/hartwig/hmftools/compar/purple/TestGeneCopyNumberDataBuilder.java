package com.hartwig.hmftools.compar.purple;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestGeneCopyNumberDataBuilder
{
    public String gene = "BRAF";
    public double minCopyNumber = 1;
    public double maxCopyNumber = 2;

    private static final Consumer<TestGeneCopyNumberDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.gene = "BRCA2";
        b.minCopyNumber = 0;
        b.maxCopyNumber = 5;
    };

    public static final TestComparableItemBuilder<TestGeneCopyNumberDataBuilder, GeneCopyNumberData> BUILDER =
            new TestComparableItemBuilder<>(
                    TestGeneCopyNumberDataBuilder::new,
                    TestGeneCopyNumberDataBuilder::build,
                    ALTERNATE_INITIALIZER
            );

    private GeneCopyNumberData build()
    {
        return new GeneCopyNumberData(
                ImmutableGeneCopyNumber.builder()
                        .geneName(gene)
                        .minCopyNumber(minCopyNumber)
                        .maxCopyNumber(maxCopyNumber)
                        .chromosome("")
                        .start(-1)
                        .end(-1)
                        .transName("")
                        .isCanonical(true)
                        .chromosomeBand("")
                        .somaticRegions(-1)
                        .minRegions(-1)
                        .minRegionStart(-1)
                        .minRegionEnd(-1)
                        .depthWindowCount(-1)
                        .gcContent(-1)
                        .minRegionStartSupport(SegmentSupport.BND)
                        .minRegionEndSupport(SegmentSupport.BND)
                        .minRegionMethod(CopyNumberMethod.UNKNOWN)
                        .minMinorAlleleCopyNumber(-1)
                        .build()
        );
    }
}
