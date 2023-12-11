package com.hartwig.hmftools.common.samtools;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_POSITION;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.N;

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

    public static int cigarBaseLength(final Cigar cigar)
    {
        return cigar.getCigarElements().stream().filter(x -> x.getOperator() != N && x.getOperator() != D).mapToInt(x -> x.getLength()).sum();
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

    public static boolean leftSoftClipped(final SAMRecord record) { return leftSoftClipped(record.getCigar()); }
    public static boolean rightSoftClipped(final SAMRecord record) { return rightSoftClipped(record.getCigar()); }

    public static boolean leftSoftClipped(final Cigar cigar)
    {
        return cigar.getCigarElements().get(0).getOperator() == CigarOperator.S;
    }

    public static boolean rightSoftClipped(final Cigar cigar)
    {
        return cigar.getCigarElements().get(cigar.getCigarElements().size() - 1).getOperator() == CigarOperator.S;
    }

    public static int leftSoftClipLength(final SAMRecord record)
    {
        CigarElement firstElement = record.getCigar().getFirstCigarElement();
        return (firstElement != null && firstElement.getOperator() == CigarOperator.S) ? firstElement.getLength() : 0;
    }

    public static int rightSoftClipLength(final SAMRecord record)
    {
        CigarElement lastElement = record.getCigar().getLastCigarElement();
        return (lastElement != null && lastElement.getOperator() == CigarOperator.S) ? lastElement.getLength() : 0;
    }

    @Nullable
    public static String leftSoftClipBases(@NotNull final SAMRecord record)
    {
        int leftClip = leftSoftClipLength(record);
        if(leftClip == 0)
            return null;

        return record.getReadString().substring(0, leftClip);
    }

    @Nullable
    public static String rightSoftClipBases(@NotNull final SAMRecord record)
    {
        int rightClip = rightSoftClipLength(record);
        if(rightClip == 0)
            return null;

        return record.getReadString().substring(record.getReadString().length() - rightClip);
    }

    public static int getUnclippedPosition(final int readStart, @NotNull final String cigarStr, final boolean forwardStrand)
    {
        return getEndPosition(readStart, cigarStr, forwardStrand, true);
    }

    public static int getEndPosition(final int readStart, @NotNull final String cigarStr, final boolean forwardStrand, boolean includeSoftClipped)
    {
        int currentPosition = readStart;
        int elementLength = 0;

        for(int i = 0; i < cigarStr.length(); ++i)
        {
            char c = cigarStr.charAt(i);
            boolean isAddItem = (c == 'D' || c == 'M' || c == 'S' || c == 'N');

            if(isAddItem)
            {
                if(forwardStrand)
                {
                    // back out the left clip if present
                    return c == 'S' ? readStart - elementLength : readStart;
                }

                if(c == 'S' && readStart == currentPosition)
                {
                    // ignore left-clip since getting the read end position
                }
                else if(c == 'S' && !includeSoftClipped)
                {
                    break;
                }
                else
                {
                    currentPosition += elementLength;
                }

                elementLength = 0;
                continue;
            }

            int digit = c - '0';
            if (digit >= 0 && digit <= 9)
            {
                elementLength = elementLength * 10 + digit;
            }
            else
            {
                elementLength = 0;
            }
        }

        // always pointing to the start of the next element, so need to move back a base
        return currentPosition - 1;
    }

    public static int getMateAlignmentEnd(final SAMRecord read)
    {
        String mateCigarStr = read.getStringAttribute(MATE_CIGAR_ATTRIBUTE);
        if(mateCigarStr == null || mateCigarStr.equals(NO_CIGAR))
            return NO_POSITION;

        return getEndPosition(read.getMateAlignmentStart(), mateCigarStr, false, false);
    }
}
