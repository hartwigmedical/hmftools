package com.hartwig.hmftools.sage.evidence;

import static java.lang.String.format;

import static com.hartwig.hmftools.sage.SageConstants.MATCHING_BASE_QUALITY;
import static com.hartwig.hmftools.sage.evidence.JitterMatch.LENGTHENED;
import static com.hartwig.hmftools.sage.evidence.JitterMatch.NONE;
import static com.hartwig.hmftools.sage.evidence.JitterMatch.SHORTENED;

import com.hartwig.hmftools.sage.common.RepeatInfo;
import com.hartwig.hmftools.sage.common.VariantReadContext;

import htsjdk.samtools.SAMRecord;

public class JitterData
{
    private int mLengthened;
    private int mShortened;

    public JitterData()
    {
        mLengthened = 0;
        mShortened = 0;
    }

    public void update(final JitterMatch jitterMatch)
    {
        if(jitterMatch == JitterMatch.LENGTHENED)
            mLengthened++;
        else if(jitterMatch == JitterMatch.SHORTENED)
            mShortened++;
    }

    public int[] summary() { return new int[] { mShortened, mLengthened, 0 }; }

    public String toString() { return format("short(%d) long(%d)", mShortened, mLengthened); }

    public static JitterMatch checkJitter(final VariantReadContext readContext, final SAMRecord record, int readVarIndex)
    {
        if(readContext.AllRepeats.isEmpty())
            return JitterMatch.NONE;

        final byte[] readBases = record.getReadBases();
        final byte[] readQuals = record.getBaseQualities();

        // try each repeat in turn, lengthening and shortening it
        for(RepeatInfo repeat : readContext.AllRepeats)
        {
            JitterMatch jitterMatch = checkRepeatJitterMatch(repeat, readContext, readVarIndex, readBases, readQuals);

            if(jitterMatch != NONE)
                return jitterMatch;
        }

        return JitterMatch.NONE;
    }

    private static JitterMatch checkRepeatJitterMatch(
            final RepeatInfo repeat, final VariantReadContext readContext, int readVarIndex, final byte[] readBases, final byte[] readQuals)
    {
        int repeatLength = repeat.repeatLength();

        for(int i = 0; i <= 1; ++i)
        {
            JitterMatch jitterType = (i == 0) ? SHORTENED : LENGTHENED;

            int flankReadIndexStart = readVarIndex - readContext.leftLength();
            int flankReadIndexEnd = readVarIndex + readContext.rightLength() - 1;

            if(jitterType == SHORTENED)
                flankReadIndexStart += repeatLength;
            else
                flankReadIndexStart -= repeatLength;

            if(flankReadIndexStart < 0 || flankReadIndexEnd >= readBases.length)
                return JitterMatch.NONE;

            int readContextIndex = 0;
            int readIndex = flankReadIndexStart;
            boolean allMatched = true;
            boolean indexAdjusted = false;

            for(; readIndex <= flankReadIndexEnd; ++readIndex, ++readContextIndex)
            {
                if(!indexAdjusted)
                {
                    if(jitterType == SHORTENED && readContextIndex == repeat.endIndex() - repeatLength)
                    {
                        indexAdjusted = true;
                        readIndex -= repeatLength;
                    }
                    else if(jitterType == LENGTHENED && readContextIndex == repeat.endIndex() + 1)
                    {
                        indexAdjusted = true;
                        readContextIndex -= repeatLength;
                    }
                }

                if(readBases[readIndex] != readContext.ReadBases[readContextIndex])
                {
                    // mismatch cannot be in the core
                    if(readContextIndex >= readContext.CoreIndexStart && readContextIndex <= readContext.CoreIndexEnd)
                    {
                        allMatched = false;
                        break;
                    }

                    // and must be low-qual
                    if(readQuals[readIndex] >= MATCHING_BASE_QUALITY)
                    {
                        allMatched = false;
                        break;
                    }
                }
            }

            if(allMatched)
                return jitterType;
        }

        return JitterMatch.NONE;
    }
}
