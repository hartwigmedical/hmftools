package com.hartwig.hmftools.common.purple;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class GeneCopyNumberTestFactory
{
    public static GeneCopyNumber createGeneCopyNumber(final String gene)
    {
        return new GeneCopyNumber(CHR_1, 0, 0, gene, Strings.EMPTY, true,
                Strings.EMPTY, 0 ,0, 0, 1, 1,
                0, 0, 0, 1.0,
                SegmentSupport.NONE, SegmentSupport.NONE, CopyNumberMethod.UNKNOWN, 0);
    }

    public static GeneCopyNumber createGeneCopyNumber(final String gene, final double minCopyNumber, final double maxCopyNumber)
    {
        return new GeneCopyNumber(CHR_1, 0, 0, gene, Strings.EMPTY, true,
                Strings.EMPTY, maxCopyNumber , minCopyNumber, 0, 1, 1,
                0, 0, 0, 1.0,
                SegmentSupport.NONE, SegmentSupport.NONE, CopyNumberMethod.UNKNOWN, minCopyNumber / 2);
    }

    public static GeneCopyNumber createGeneCopyNumber(
            final String chromosome, final String gene, final double minCopyNumber, final double maxCopyNumber)
    {
        return new GeneCopyNumber(chromosome, 0, 0, gene, Strings.EMPTY, true,
                Strings.EMPTY, maxCopyNumber , minCopyNumber, 0, 1, 1,
                0, 0, 0, 1.0,
                SegmentSupport.NONE, SegmentSupport.NONE, CopyNumberMethod.UNKNOWN, minCopyNumber / 2);
    }

    /*
    @NotNull
    public static ImmutableGeneCopyNumber.Builder builder()
    {
        return ImmutableGeneCopyNumber.builder()
                .start(0)
                .end(0)
                .geneName(Strings.EMPTY)
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .minRegionStart(0)
                .minRegionStartSupport(SegmentSupport.NONE)
                .minRegionEnd(0)
                .minRegionEndSupport(SegmentSupport.NONE)
                .minRegionMethod(CopyNumberMethod.UNKNOWN)
                .minRegions(1)
                .somaticRegions(1)
                .minCopyNumber(0D)
                .maxCopyNumber(0D)
                .transName(Strings.EMPTY)
                .isCanonical(true)
                .minMinorAlleleCopyNumber(0)
                .depthWindowCount(0)
                .gcContent(1.0);
    }
    */
}
