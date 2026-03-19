package com.hartwig.hmftools.common.purple;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class GeneCopyNumberTestFactory
{
    public static GeneCopyNumber createGeneCopyNumber(final String gene, final double minCopyNumber, final double maxCopyNumber)
    {
        return new GeneCopyNumber(CHR_1, 0, 0, gene, Strings.EMPTY, true,
                "p36.13", maxCopyNumber , minCopyNumber, 0, 1, 1,
                0, 0, 0, 1.0,
                SegmentSupport.NONE, SegmentSupport.NONE, CopyNumberMethod.UNKNOWN, minCopyNumber / 2);
    }

    public static GeneCopyNumber createGeneCopyNumber(
            final String chromosome, final String gene, final double minCopyNumber, final double maxCopyNumber)
    {
        return new GeneCopyNumber(chromosome, 0, 0, gene, Strings.EMPTY, true,
                "p36.13", maxCopyNumber , minCopyNumber, 0, 1, 1,
                0, 0, 0, 1.0,
                SegmentSupport.NONE, SegmentSupport.NONE, CopyNumberMethod.UNKNOWN, minCopyNumber / 2);
    }
}
