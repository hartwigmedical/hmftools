package com.hartwig.hmftools.esvee.read;

import static com.hartwig.hmftools.esvee.SvConstants.INDEL_TO_SC_MAX_EDGE_DISTANCE;
import static com.hartwig.hmftools.esvee.SvConstants.INDEL_TO_SC_MIN_SIZE_SOFTCLIP;
import static com.hartwig.hmftools.esvee.SvConstants.POLY_G_TRIM_LENGTH;
import static com.hartwig.hmftools.esvee.common.BaseType.G;
import static com.hartwig.hmftools.esvee.common.BaseType.C;

import static htsjdk.samtools.CigarOperator.I;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class ReadAdjustments
{

    public ReadAdjustments()
    {

    }

    public static boolean convertEdgeIndelsToSoftClip(final Read read)
    {
        if(read.cigarElements().size() < 3)
            return false;

        int leftSoftClipLength = calcIndelToSoftClipLength(
                read.cigarElements().get(0), read.cigarElements().get(1), read.cigarElements().get(2));

        int rightSoftClipLength = 0;

        if(leftSoftClipLength == 0 || read.cigarElements().size() >= 5) // cannot adjust the same I from each side so re-check length
        {
            int lastIndex = read.cigarElements().size() - 1;

            rightSoftClipLength = calcIndelToSoftClipLength(
                    read.cigarElements().get(lastIndex), read.cigarElements().get(lastIndex - 1), read.cigarElements().get(lastIndex - 2));
        }

        if(leftSoftClipLength > 0 || rightSoftClipLength > 0)
        {
            read.convertEdgeIndelToSoftClip(leftSoftClipLength, rightSoftClipLength);
            return true;
        }

        return false;
    }

    private static int calcIndelToSoftClipLength(final CigarElement edge, final CigarElement inside, final CigarElement next)
    {
        if(edge.getOperator() != CigarOperator.M || edge.getLength() >= INDEL_TO_SC_MAX_EDGE_DISTANCE)
            return 0;

        if(!inside.getOperator().isIndel())
            return 0;

        if(next.getOperator() != CigarOperator.M)
            return 0;

        if(inside.getLength() < INDEL_TO_SC_MIN_SIZE_SOFTCLIP)
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


}
