package com.hartwig.hmftools.compar.purple;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestCopyNumberDataBuilder
{
    public String chromosome = "chr1";
    public int positionStart = 1;
    public int positionEnd = 100000000;
    public double copyNumber = 2;
    public double majorAlleleCopyNumber = 1;
    public CopyNumberMethod method = CopyNumberMethod.BAF_WEIGHTED;
    public String comparisonChromosomeStart = "chr1";
    public String comparisonChromosomeEnd = "chr1";
    public int comparisonPositionStart = 1;
    public int comparisonPositionEnd = 100000000;

    private static final Consumer<TestCopyNumberDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.chromosome = "chr8";
        b.positionStart = 10000;
        b.positionEnd = 30000;
        b.copyNumber = 3;
        b.majorAlleleCopyNumber = 3;
        b.method = CopyNumberMethod.STRUCTURAL_VARIANT;
        b.comparisonChromosomeStart = "chr8";
        b.comparisonChromosomeEnd = "chr8";
        b.comparisonPositionStart = 10000;
        b.comparisonPositionEnd = 30000;
    };

    public static final TestComparableItemBuilder<TestCopyNumberDataBuilder, CopyNumberData> BUILDER =
            new TestComparableItemBuilder<>(TestCopyNumberDataBuilder::new, TestCopyNumberDataBuilder::build, ALTERNATE_INITIALIZER);

    private CopyNumberData build()
    {
        return new CopyNumberData(
                chromosome,
                positionStart,
                positionEnd,
                copyNumber,
                majorAlleleCopyNumber,
                method,
                new BasePosition(comparisonChromosomeStart, comparisonPositionStart),
                new BasePosition(comparisonChromosomeEnd, comparisonPositionEnd)
        );
    }
}
