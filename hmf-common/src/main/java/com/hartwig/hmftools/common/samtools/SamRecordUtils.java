package com.hartwig.hmftools.common.samtools;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public final class SamRecordUtils
{
    public static final int PHRED_OFFSET = 33;

    public static int leftSoftClip(@NotNull final SAMRecord record)
    {
        Cigar cigar = record.getCigar();
        if(cigar.numCigarElements() > 0)
        {
            CigarElement firstElement = cigar.getCigarElement(0);
            if(firstElement.getOperator() == CigarOperator.S)
            {
                return firstElement.getLength();
            }
        }

        return 0;
    }

    public static int rightSoftClip(@NotNull final SAMRecord record)
    {
        Cigar cigar = record.getCigar();
        if(cigar.numCigarElements() > 0)
        {
            CigarElement lastElement = cigar.getCigarElement(cigar.numCigarElements() - 1);
            if(lastElement.getOperator() == CigarOperator.S)
            {
                return lastElement.getLength();
            }
        }

        return 0;
    }

    @Nullable
    public static String leftSoftClipBases(@NotNull final SAMRecord record)
    {
        int leftClip = leftSoftClip(record);
        if (leftClip == 0)
        {
            return null;
        }
        return record.getReadString().substring(0, leftClip);
    }

    @Nullable
    public static String rightSoftClipBases(@NotNull final SAMRecord record)
    {
        int rightClip = rightSoftClip(record);
        if (rightClip == 0)
        {
            return null;
        }
        return record.getReadString().substring(record.getReadString().length() - rightClip);
    }

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
