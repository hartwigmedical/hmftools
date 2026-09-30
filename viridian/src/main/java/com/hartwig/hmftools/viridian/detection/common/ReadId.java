package com.hartwig.hmftools.viridian.detection.common;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

// Uniquely identifies a read.
public record ReadId(
        String name,
        // 1 or 2, or 0 where the read is unpaired
        int mateIndex
)
{
    private static final int UNPAIRED = 0;
    private static final char SUFFIX_SEPARATOR = '/';
    private static final int SUFFIX_LENGTH = 2;

    public ReadId
    {
        if(name.isEmpty())
        {
            throw new IllegalArgumentException("Invalid name");
        }
        if(mateIndex < UNPAIRED || mateIndex > 2)
        {
            throw new IllegalArgumentException("Invalid mateIndex: " + mateIndex);
        }
    }

    public static ReadId from(SAMRecord record)
    {
        int mateIndex = UNPAIRED;
        if(record.getReadPairedFlag())
        {
            mateIndex = record.getFirstOfPairFlag() ? 1 : 2;
        }
        return new ReadId(record.getReadName(), mateIndex);
    }

    public static ReadId parse(String string)
    {
        int suffixStart = string.length() - SUFFIX_LENGTH;
        if(suffixStart > 0 && string.charAt(suffixStart) == SUFFIX_SEPARATOR)
        {
            int mateIndex = string.charAt(string.length() - 1) - '0';
            if(mateIndex == 1 || mateIndex == 2)
            {
                return new ReadId(string.substring(0, suffixStart), mateIndex);
            }
        }
        return new ReadId(string, UNPAIRED);
    }

    @NotNull
    @Override
    public String toString()
    {
        return mateIndex == UNPAIRED ? name : name + SUFFIX_SEPARATOR + mateIndex;
    }
}
