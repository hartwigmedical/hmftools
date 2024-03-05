package com.hartwig.hmftools.esvee.read;

import static com.hartwig.hmftools.esvee.SvConstants.INDEL_TO_SC_MAX_EDGE_DISTANCE;
import static com.hartwig.hmftools.esvee.SvConstants.INDEL_TO_SC_MAX_SIZE_SOFTCLIP;
import static com.hartwig.hmftools.esvee.SvConstants.INDEL_TO_SC_MIN_SIZE_SOFTCLIP;
import static com.hartwig.hmftools.esvee.SvConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.SvConstants.LOW_BASE_TRIM_PERC;
import static com.hartwig.hmftools.esvee.SvConstants.POLY_G_TRIM_LENGTH;
import static com.hartwig.hmftools.esvee.common.BaseType.G;
import static com.hartwig.hmftools.esvee.common.BaseType.C;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public final class ReadAdjustments
{
    public static boolean convertEdgeIndelsToSoftClip(final Read read)
    {
        return convertEdgeIndelsToSoftClip(read, INDEL_TO_SC_MAX_SIZE_SOFTCLIP, INDEL_TO_SC_MAX_EDGE_DISTANCE);
    }

    public static boolean convertEdgeIndelsToSoftClip(final Read read, final int maxIndelLength, final int maxReadEdgeDistance)
    {
        if(read.cigarElements().size() < 3)
            return false;

        int leftSoftClipLength = calcIndelToSoftClipLength(
                read.cigarElements().get(0), read.cigarElements().get(1), read.cigarElements().get(2), maxIndelLength, maxReadEdgeDistance);

        int rightSoftClipLength = 0;

        if(leftSoftClipLength == 0 || read.cigarElements().size() >= 5) // cannot adjust the same I from each side so re-check length
        {
            int lastIndex = read.cigarElements().size() - 1;

            rightSoftClipLength = calcIndelToSoftClipLength(
                    read.cigarElements().get(lastIndex), read.cigarElements().get(lastIndex - 1), read.cigarElements().get(lastIndex - 2),
                    maxIndelLength, maxReadEdgeDistance);
        }

        if(leftSoftClipLength > 0 || rightSoftClipLength > 0)
        {
            read.convertEdgeIndelToSoftClip(leftSoftClipLength, rightSoftClipLength);
            return true;
        }

        return false;
    }

    private static int calcIndelToSoftClipLength(
            final CigarElement edge, final CigarElement inside, final CigarElement next,
            final int maxIndelLength, final int maxReadEdgeDistance)
    {
        if(edge.getOperator() != CigarOperator.M || edge.getLength() >= maxReadEdgeDistance)
            return 0;

        if(!inside.getOperator().isIndel())
            return 0;

        if(next.getOperator() != CigarOperator.M)
            return 0;

        if(inside.getLength() < INDEL_TO_SC_MIN_SIZE_SOFTCLIP || inside.getLength() > maxIndelLength)
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
