package com.hartwig.hmftools.esvee.read;

import static com.hartwig.hmftools.esvee.AssemblyConstants.INDEL_TO_SC_MAX_SIZE_SOFTCLIP;
import static com.hartwig.hmftools.esvee.AssemblyConstants.INDEL_TO_SC_MIN_SIZE_SOFTCLIP;
import static com.hartwig.hmftools.esvee.AssemblyConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.AssemblyConstants.LOW_BASE_TRIM_PERC;
import static com.hartwig.hmftools.esvee.AssemblyConstants.POLY_G_TRIM_LENGTH;
import static com.hartwig.hmftools.esvee.types.BaseType.G;
import static com.hartwig.hmftools.esvee.types.BaseType.C;

import static htsjdk.samtools.CigarOperator.M;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public final class ReadAdjustments
{
    public static boolean convertEdgeIndelsToSoftClip(final Read read)
    {
        return convertEdgeIndelsToSoftClip(read, INDEL_TO_SC_MIN_SIZE_SOFTCLIP, INDEL_TO_SC_MAX_SIZE_SOFTCLIP, true);
    }

    public static boolean convertEdgeIndelsToSoftClip(
            final Read read, final int minIndelLength, final int maxIndelLength, boolean allowDoubleConversion)
    {
        if(read.cigarElements().size() < 3)
            return false;

        int leftSoftClipLength = calcIndelToSoftClipLength(
                read.cigarElements().get(0), read.cigarElements().get(1), read.cigarElements().get(2),
                minIndelLength, maxIndelLength);

        int leftEdgeDistance = leftSoftClipLength > 0 ? read.cigarElements().get(0).getLength() : 0;

        int lastIndex = read.cigarElements().size() - 1;

        int rightSoftClipLength = calcIndelToSoftClipLength(
                read.cigarElements().get(lastIndex), read.cigarElements().get(lastIndex - 1), read.cigarElements().get(lastIndex - 2),
                minIndelLength, maxIndelLength);

        int rightEdgeDistance = rightSoftClipLength > 0 ? read.cigarElements().get(lastIndex).getLength() : 0;

        if(leftSoftClipLength > 0 && rightSoftClipLength > 0 && !allowDoubleConversion)
        {
            if(leftEdgeDistance < rightEdgeDistance)
                rightSoftClipLength = 0;
            else
                leftSoftClipLength = 0;
        }

        if(leftSoftClipLength > 0 || rightSoftClipLength > 0)
        {
            if(allowDoubleConversion)
                read.setIndelUnclippedBounds(leftSoftClipLength, rightSoftClipLength);
            else
                read.convertEdgeIndelToSoftClip(leftSoftClipLength, rightSoftClipLength);

            return true;
        }

        return false;
    }

    public static void convertIndelSoftClip(final Read read, boolean onLeft)
    {
        if(read.cigarElements().size() < 3)
            return;

        if(onLeft)
        {
            int leftSoftClipLength = calcIndelToSoftClipLength(
                    read.cigarElements().get(0), read.cigarElements().get(1), read.cigarElements().get(2),
                    INDEL_TO_SC_MIN_SIZE_SOFTCLIP, INDEL_TO_SC_MAX_SIZE_SOFTCLIP);

            read.convertEdgeIndelToSoftClip(leftSoftClipLength, 0);
        }
        else
        {
            int lastIndex = read.cigarElements().size() - 1;

            int rightSoftClipLength = calcIndelToSoftClipLength(
                    read.cigarElements().get(lastIndex), read.cigarElements().get(lastIndex - 1), read.cigarElements().get(lastIndex - 2),
                    INDEL_TO_SC_MIN_SIZE_SOFTCLIP, INDEL_TO_SC_MAX_SIZE_SOFTCLIP);

            read.convertEdgeIndelToSoftClip(0, rightSoftClipLength);
        }
    }

    private static int calcIndelToSoftClipLength(
            final CigarElement edge, final CigarElement inside, final CigarElement next,
            final int minIndelLength, final int maxIndelLength)
    {
        if(edge.getOperator() != M)
            return 0;

        if(!inside.getOperator().isIndel())
            return 0;

        if(next.getOperator() != M)
            return 0;

        if(inside.getLength() < minIndelLength || inside.getLength() > maxIndelLength)
            return 0;

        return inside.getOperator() != CigarOperator.D ? edge.getLength() + inside.getLength() : edge.getLength();
    }

    public static boolean trimPolyGSequences(final Read read)
    {
        if(read.isUnmapped() || !read.isPairedRead())
            return false;

        int trailingGCount = 0;

        if(read.positiveStrand())
        {
            for(int i = read.basesLength() - 1; i >= 0; --i)
            {
                if(read.getBases()[i] == G.Byte)
                    trailingGCount++;
                else
                    break;
            }
        }
        else
        {
            for(int i = 0; i < read.basesLength(); ++i)
            {
                if(read.getBases()[i] == C.Byte)
                    trailingGCount++;
                else
                    break;
            }
        }

        if(trailingGCount < POLY_G_TRIM_LENGTH)
            return false;

        read.trimBases(trailingGCount, read.negativeStrand());

        return true;
    }

    public static boolean trimLowQualBases(final Read read)
    {
        boolean fromStart = read.negativeStrand();
        int scBaseCount = fromStart ? read.leftClipLength() : read.rightClipLength();

        if(scBaseCount == 0)
            return false;

        int baseIndex = fromStart ? 0 : read.basesLength() - 1;

        int lowQualCount = 0;
        int lastLowQualPercIndex = 0;

        for(int i = 1; i <= scBaseCount; ++i)
        {
            if(read.getBaseQuality()[baseIndex] < LOW_BASE_QUAL_THRESHOLD)
            {
                lowQualCount++;

                if(lowQualCount / (double)i >= LOW_BASE_TRIM_PERC)
                    lastLowQualPercIndex = i;
            }

            if(fromStart)
                ++baseIndex;
            else
                --baseIndex;
        }

        if(lastLowQualPercIndex == 0)
            return false;

        read.trimBases(lastLowQualPercIndex, read.negativeStrand());

        return true;
    }
}
