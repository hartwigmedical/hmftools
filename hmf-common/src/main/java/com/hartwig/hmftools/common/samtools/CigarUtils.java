package com.hartwig.hmftools.common.samtools;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public final class CigarUtils
{
    public static Cigar cigarFromStr(final String cigarStr)
    {
        Cigar cigar = new Cigar();

        int length = 0;
        for (int i = 0; i < cigarStr.length(); ++i)
        {
            char c = cigarStr.charAt(i);
            int digit = c - '0';

            if (digit >= 0 && digit <= 9)
            {
                length = length * 10 + digit;
            }
            else
            {
                CigarOperator operator = CigarOperator.characterToEnum(c);
                cigar.add(new CigarElement(length, operator));
                length = 0;
            }
        }

        return cigar;
    }

    public static int calcCigarLength(final String cigarStr)
    {
        int baseLength = 0;
        int currentElementLength = 0;
        for (int i = 0; i < cigarStr.length(); ++i)
        {
            char c = cigarStr.charAt(i);
            boolean isAddItem = (c == 'D' || c == 'M');

            if(isAddItem)
            {
                baseLength += currentElementLength;
                currentElementLength = 0;
                continue;
            }
            int digit = c - '0';
            if (digit >= 0 && digit <= 9)
            {
                currentElementLength = currentElementLength * 10 + digit;
            }
            else
            {
                // we are in one of the ignored items such as S, H, N etc
                currentElementLength = 0;
            }
        }

        return baseLength;
    }

    public static int leftSoftClip(@NotNull final SAMRecord record)
    {
        CigarElement firstElement = record.getCigar().getFirstCigarElement();
        return (firstElement != null && firstElement.getOperator() == CigarOperator.S) ? firstElement.getLength() : 0;
    }

    public static int rightSoftClip(@NotNull final SAMRecord record)
    {
        CigarElement lastElement = record.getCigar().getLastCigarElement();
        return (lastElement != null && lastElement.getOperator() == CigarOperator.S) ? lastElement.getLength() : 0;
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
}
