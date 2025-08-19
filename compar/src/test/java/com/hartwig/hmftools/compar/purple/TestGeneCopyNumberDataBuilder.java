package com.hartwig.hmftools.compar.purple;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

import org.apache.logging.log4j.util.Strings;

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
        GeneCopyNumber geneCopyNumber = new GeneCopyNumber("", 0, 0, gene, Strings.EMPTY, true,
                    Strings.EMPTY, maxCopyNumber , minCopyNumber, 0, 1, 1,
                    0, 0, 0, 1.0,
                    SegmentSupport.BND, SegmentSupport.BND, CopyNumberMethod.UNKNOWN);

        return new GeneCopyNumberData(geneCopyNumber);
    }
}
