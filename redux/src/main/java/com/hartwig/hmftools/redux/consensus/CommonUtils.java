package com.hartwig.hmftools.redux.consensus;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASE_BYTES;
import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_N_BYTE;
import static com.hartwig.hmftools.redux.consensus.BaseQualPair.NO_BASE;

import htsjdk.samtools.SAMRecord;

public final class CommonUtils
{
    public static boolean hasValidBases(final SAMRecord record)
    {
        for(int i = 0; i < record.getReadBases().length; ++i)
        {
            if(record.getReadBases()[i] == NO_BASE)
                return false;
        }

        return true;
    }

    public static int findMostCommonBaseCount(final int[] baseCounts)
    {
        int maxCount = 0;

        for(int b = 0; b < DNA_BASE_BYTES.length; ++b)
        {
            maxCount = max(maxCount, baseCounts[b]);
        }

        return maxCount;
    }

    public static byte findMostCommonBase(final int[] baseCounts, final byte refBase, final int maxCount)
    {
        byte maxBase = DNA_N_BYTE;

        for(int b = 0; b < DNA_BASE_BYTES.length; ++b)
        {
            if(baseCounts[b] == maxCount)
            {
                if(DNA_BASE_BYTES[b] == refBase || maxBase == DNA_N_BYTE)
                    maxBase = DNA_BASE_BYTES[b];
            }
        }

        return maxBase;
    }
}
