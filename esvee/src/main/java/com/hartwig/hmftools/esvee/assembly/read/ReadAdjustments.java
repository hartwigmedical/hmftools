package com.hartwig.hmftools.esvee.assembly.read;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.sv.LineElements.LINE_POLY_AT_REQ;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_POLY_AT_TEST_LEN;
import static com.hartwig.hmftools.esvee.AssemblyConstants.INDEL_TO_SC_MAX_SIZE_SOFTCLIP;
import static com.hartwig.hmftools.esvee.AssemblyConstants.INDEL_TO_SC_MIN_SIZE_SOFTCLIP;
import static com.hartwig.hmftools.esvee.AssemblyConstants.LOW_BASE_TRIM_PERC;
import static com.hartwig.hmftools.esvee.AssemblyConstants.POLY_G_TRIM_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.findLineSequenceBase;
import static com.hartwig.hmftools.esvee.common.SvConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.assembly.types.BaseType.G;
import static com.hartwig.hmftools.esvee.assembly.types.BaseType.C;

import static htsjdk.samtools.CigarOperator.M;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public final class ReadAdjustments
{
    public static boolean convertEdgeIndelsToSoftClip(final Read read)
    {
        return convertEdgeIndelsToSoftClip(read, INDEL_TO_SC_MIN_SIZE_SOFTCLIP, INDEL_TO_SC_MAX_SIZE_SOFTCLIP);
    }

    public static boolean convertEdgeIndelsToSoftClip(final Read read, final int minIndelLength, final int maxIndelLength)
    {
        if(read.cigarElements().size() < 3)
            return false;

        int leftSoftClipLength = calcIndelToSoftClipLength(
                read.cigarElements().get(0), read.cigarElements().get(1), read.cigarElements().get(2),
                minIndelLength, maxIndelLength);

        int lastIndex = read.cigarElements().size() - 1;

        int rightSoftClipLength = calcIndelToSoftClipLength(
                read.cigarElements().get(lastIndex), read.cigarElements().get(lastIndex - 1), read.cigarElements().get(lastIndex - 2),
                minIndelLength, maxIndelLength);

        if(leftSoftClipLength > 0 || rightSoftClipLength > 0)
        {
            read.setIndelUnclippedBounds(leftSoftClipLength, rightSoftClipLength);
            return true;
        }

        return false;
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

    public static boolean trimLowQualSoftClipBases(final Read read)
    {
        boolean fromStart = read.negativeStrand();
        int scBaseCount = fromStart ? read.leftClipLength() : read.rightClipLength();

        if(scBaseCount == 0)
            return false;

        // first establish the 5' end of the soft-clip satisfies LINE criteria
        int lineExclusionLength = 0;

        if(scBaseCount >= LINE_POLY_AT_REQ)
        {
            int scIndexStart, scIndexEnd;
            int lineTestLength = min(scBaseCount, LINE_POLY_AT_TEST_LEN);

            if(fromStart)
            {
                scIndexEnd = scBaseCount - 1;
                scIndexStart = scIndexEnd - lineTestLength + 1;
            }
            else
            {
                scIndexStart = read.basesLength() - scBaseCount;
                scIndexEnd = scIndexStart + lineTestLength - 1;
            }

            Byte lineBase = findLineSequenceBase(read.getBases(), scIndexStart, scIndexEnd);

            if(lineBase != null)
            {
                read.markLineTail();

                lineExclusionLength = lineTestLength;

                // find the most 3' index for the observed line base

                int baseIndex = fromStart ? scIndexStart : scIndexEnd;
                int baseCheck = lineTestLength - LINE_POLY_AT_REQ;

                while(baseCheck > 0)
                {
                    if(read.getBases()[baseIndex] != lineBase)
                        --lineExclusionLength;

                    if(fromStart)
                        ++baseIndex;
                    else
                        --baseIndex;

                    --baseCheck;
                }
            }
        }

        int baseIndex = fromStart ? 0 : read.basesLength() - 1;

        int lowQualCount = 0;
        int lastLowQualPercIndex = 0;

        for(int i = 1; i <= scBaseCount - lineExclusionLength; ++i)
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

        read.trimBases(lastLowQualPercIndex, fromStart);

        return true;
    }

    public static boolean trimLowQualBases(final Read read)
    {
        boolean fromStart = read.negativeStrand();

        int baseIndex = fromStart ? 0 : read.basesLength() - 1;

        int lowQualCount = 0;

        while(baseIndex != read.basesLength() / 2) // at most half the read
        {
            if(read.getBaseQuality()[baseIndex] < LOW_BASE_QUAL_THRESHOLD)
                lowQualCount++;
            else
                break;

            if(fromStart)
                ++baseIndex;
            else
                --baseIndex;
        }

        if(lowQualCount == 0)
            return false;

        read.trimBases(lowQualCount, fromStart);

        return true;
    }
}
