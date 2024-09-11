package com.hartwig.hmftools.esvee.assembly.read;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.sv.LineElements.LINE_BASE_A;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_BASE_T;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_POLY_AT_REQ;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_POLY_AT_TEST_LEN;
import static com.hartwig.hmftools.esvee.AssemblyConstants.INDEL_TO_SC_MAX_SIZE_SOFTCLIP;
import static com.hartwig.hmftools.esvee.AssemblyConstants.INDEL_TO_SC_MIN_SIZE_SOFTCLIP;
import static com.hartwig.hmftools.esvee.AssemblyConstants.LOW_BASE_TRIM_PERC;
import static com.hartwig.hmftools.esvee.AssemblyConstants.POLY_G_TRIM_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.LineUtils.findBaseRepeatCount;
import static com.hartwig.hmftools.esvee.assembly.LineUtils.findLineSequenceCount;
import static com.hartwig.hmftools.esvee.common.CommonUtils.belowMinQual;
import static com.hartwig.hmftools.esvee.common.SvConstants.LINE_REF_BASE_REPEAT_FACTOR;
import static com.hartwig.hmftools.esvee.assembly.types.BaseType.G;
import static com.hartwig.hmftools.esvee.assembly.types.BaseType.C;

import static htsjdk.samtools.CigarOperator.M;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

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

    public static void markLineSoftClips(final Read read)
    {
        for(int i = 0; i <= 1; ++i)
        {
            boolean fromStart = (i == 0);
            int scBaseCount = fromStart ? read.leftClipLength() : read.rightClipLength();

            if(scBaseCount == 0)
                continue;

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

                byte lineBase = fromStart ? LINE_BASE_A : LINE_BASE_T;
                int lineBaseCount = findLineSequenceCount(read.getBases(), scIndexStart, scIndexEnd, lineBase);

                if(lineBaseCount == 0)
                    continue;

                // test that the LINE sequence doesn't extend a long repeat of the same base
                int refBaseIndex = fromStart ? scBaseCount : read.basesLength() - scBaseCount - 1;

                int refRepeatLength = findBaseRepeatCount(read.getBases(), refBaseIndex, fromStart, lineBase);

                if(lineBaseCount >= refRepeatLength * LINE_REF_BASE_REPEAT_FACTOR)
                {
                    read.markLineTail();
                    return;
                }
            }
        }
    }

    public static boolean filterLowQualRead(final SAMRecord read)
    {
        int baseLength = read.getReadBases().length;
        int qualCountThreshold = baseLength / 2 + 1;
        int lowQualCount = 0;

        for(int i = 0; i < baseLength; ++i)
        {
            if(belowMinQual(read.getBaseQualities()[i]))
            {
                ++lowQualCount;

                if(lowQualCount >= qualCountThreshold)
                    return true;
            }
            else
            {
                // exit early if majority will be high-qual
                int highQualCount = i + 1 - lowQualCount;
                if(highQualCount >= qualCountThreshold)
                    return false;
            }
        }

        return false;
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

            byte lineBase = fromStart ? LINE_BASE_A : LINE_BASE_T;

            if(read.hasLineTail())
            {
                lineExclusionLength = lineTestLength;

                // find the outermost index for the observed line base
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
            if(belowMinQual(read.getBaseQuality()[baseIndex]))
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

    public synchronized static void trimLowQualBases(final Read read)
    {
        if(read.lowQualTrimmed())
            return;

        boolean fromStart = read.negativeStrand();

        int baseLength = read.basesLength();
        int baseIndex = fromStart ? 0 : baseLength - 1;

        int lowQualCount = 0;
        int lastLowQualPercIndex = 0;
        int checkedBases = 0;

        while(baseIndex >= 0 && baseIndex < baseLength)
        {
            ++checkedBases;

            if(belowMinQual(read.getBaseQuality()[baseIndex]))
            {
                lowQualCount++;

                if(lowQualCount / (double)checkedBases >= LOW_BASE_TRIM_PERC)
                    lastLowQualPercIndex = checkedBases;
            }

            if(fromStart)
                ++baseIndex;
            else
                --baseIndex;
        }

        if(lastLowQualPercIndex > 0)
        {
            read.trimBases(lastLowQualPercIndex, fromStart);
            read.markLowQualTrimmed();
        }
    }

    public synchronized static void trimLowQualBasesStrict(final Read read)
    {
        // breaks at the first non-low-qual base from the 3' end
        if(read.lowQualTrimmed())
            return;

        boolean fromStart = read.negativeStrand();

        int baseLength = read.basesLength();
        int baseIndex = fromStart ? 0 : baseLength - 1;

        int lowQualCount = 0;

        while(baseIndex >= 0 && baseIndex < baseLength)
        {
            if(belowMinQual(read.getBaseQuality()[baseIndex]))
            {
                lowQualCount++;
                break;
            }

            if(fromStart)
                ++baseIndex;
            else
                --baseIndex;
        }

        if(lowQualCount > 0)
        {
            read.trimBases(lowQualCount, fromStart);
            read.markLowQualTrimmed();
        }
    }
}
