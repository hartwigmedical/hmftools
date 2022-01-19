package com.hartwig.hmftools.common.samtools;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public final class SamRecordUtils
{
    public static final int PHRED_OFFSET = 33;

    public static int getBaseQuality(final char quality)
    {
        return quality - PHRED_OFFSET;
    }

    public static int getBaseQuality(@NotNull final SAMRecord record, int readPosition)
    {
        return getAvgBaseQuality(record, readPosition, 1);
    }

    public static int getAvgBaseQuality(@NotNull final SAMRecord record, int readPosition, int length)
    {
        assert (readPosition > 0);

        int score = 0;
        final String baseQualities = record.getBaseQualityString();
        for(int index = readPosition - 1; index < Math.min(readPosition - 1 + length, baseQualities.length()); index++)
        {
            int baseScore = getBaseQuality(baseQualities.charAt(index));
            score += baseScore;
        }
        return score / length;
    }
}
