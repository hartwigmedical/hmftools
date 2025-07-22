package com.hartwig.hmftools.esvee.assembly.read;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.sv.LineElements.LINE_BASE_A;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_BASE_T;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_POLY_AT_REQ;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_POLY_AT_TEST_LEN;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.LOW_BASE_TRIM_PERC;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.POLY_G_TRIM_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.LineUtils.findBaseRepeatCount;
import static com.hartwig.hmftools.esvee.assembly.LineUtils.findLineSequenceCount;
import static com.hartwig.hmftools.esvee.common.CommonUtils.belowMinQual;
import static com.hartwig.hmftools.esvee.common.SvConstants.LINE_REF_BASE_REPEAT_FACTOR;
import static com.hartwig.hmftools.esvee.assembly.types.BaseType.G;
import static com.hartwig.hmftools.esvee.assembly.types.BaseType.C;

public final class ReadAdjustments
{
    public static boolean trimPolyGSequences(final Read read)
    {
        // poly-G may be trimmed from the 3' read end
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

        int readIndexStart, readIndexEnd;

        int scCheckLength = scBaseCount - lineExclusionLength;

        if(fromStart)
        {
            readIndexStart = 0;
            readIndexEnd = readIndexStart + scCheckLength - 1;
        }
        else
        {
            readIndexEnd = read.basesLength() - 1;
            readIndexStart = readIndexEnd - scCheckLength + 1;

        }

        int trimCount = findLowBaseQualTrimCount(read, readIndexStart, readIndexEnd);

        if(trimCount <= 0)
            return false;

        read.trimBases(trimCount, fromStart);
        return true;
    }

    public static int findLowBaseQualTrimCount(final Read read, int readIndexStart, int readIndexEnd)
    {
        boolean fromStart = read.negativeStrand();

        double lowestScore = 0;
        double currentScore = 0;
        int lastLowestScoreIndex = -1;

        int baseIndex = fromStart ? readIndexStart : readIndexEnd;

        while(baseIndex >= readIndexStart && baseIndex <= readIndexEnd)
        {
            if(belowMinQual(read.getBaseQuality()[baseIndex]))
            {
                currentScore -= LOW_QUAL_SCORE;

                if(currentScore <= lowestScore)
                {
                    lastLowestScoreIndex = baseIndex;
                    lowestScore = currentScore;
                }
            }
            else
            {
                ++currentScore;
            }

            if(fromStart)
                ++baseIndex;
            else
                --baseIndex;
        }

        if(lastLowestScoreIndex < 0)
            return 0;

        return fromStart ? lastLowestScoreIndex - readIndexStart + 1 : readIndexEnd - lastLowestScoreIndex + 1;
    }

    protected static final double LOW_QUAL_SCORE = 1 / LOW_BASE_TRIM_PERC - 1;
}
