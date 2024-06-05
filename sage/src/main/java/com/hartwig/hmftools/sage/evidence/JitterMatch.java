package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.sage.SageConstants.MATCHING_BASE_QUALITY;

import com.hartwig.hmftools.sage.common.RepeatInfo;
import com.hartwig.hmftools.sage.common.VariantReadContext;

import htsjdk.samtools.SAMRecord;

public enum JitterMatch
{
    SHORTENED,
    LENGTHENED,
    NONE;

    public static JitterMatch checkJitter(final VariantReadContext readContext, final SAMRecord record, int readVarIndex)
    {
        if(readContext.AllRepeats.isEmpty())
            return JitterMatch.NONE;

        final byte[] readBases = record.getReadBases();
        final byte[] readQuals = record.getBaseQualities();

        // try each repeat covering the read context in turn
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

        int readVarIndexOffset = readVarIndex - readContext.VarIndex;;
        int relativeRepeatIndex = repeat.Index + readVarIndexOffset;
        int relativeRepeatIndexEnd = repeat.endIndex() + readVarIndexOffset;

        // cannot test a repeat if the read does not fully cover it
        if(relativeRepeatIndex < 0 || relativeRepeatIndexEnd >= readBases.length)
            return JitterMatch.NONE;

        // for each read try shortening and then lengthening it
        for(int i = 0; i <= 1; ++i)
        {
            JitterMatch jitterType = (i == 0) ? SHORTENED : LENGTHENED;

            int flankReadIndexStart = readVarIndex - readContext.leftLength();
            int flankReadIndexEnd = readVarIndex + readContext.rightLength() - 1;

            if(relativeRepeatIndex < readVarIndex)
            {
                // factor in repeats before the variant read index by shifting the implied flank start
                if(jitterType == SHORTENED)
                    flankReadIndexStart += repeatLength;
                else
                    flankReadIndexStart -= repeatLength;
            }

            // at least one flank must be fully present and the core must be covered
            int coreReadIndexStart = flankReadIndexStart + readContext.leftFlankLength();
            int coreReadIndexEnd = flankReadIndexEnd - readContext.rightFlankLength();

            if(coreReadIndexStart < 0 || coreReadIndexEnd >= readBases.length)
                return JitterMatch.NONE;

            int readContextIndex = 0;
            int readIndex = flankReadIndexStart;
            boolean allMatched = true;
            boolean indexAdjusted = false;

            for(; readIndex <= flankReadIndexEnd; ++readIndex, ++readContextIndex)
            {
                if(readIndex < 0)
                    continue;

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

                if(readIndex >= readBases.length || readContextIndex >= readContext.ReadBases.length)
                    break;

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
