package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.sage.SageConstants.DOUBLE_JITTER_REPEAT_COUNT;
import static com.hartwig.hmftools.sage.SageConstants.MATCHING_BASE_QUALITY;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.sage.common.ReadContextMatcher;
import com.hartwig.hmftools.sage.common.RepeatInfo;
import com.hartwig.hmftools.sage.common.VariantReadContext;

import htsjdk.samtools.SAMRecord;

public enum JitterMatch
{
    SHORTENED,
    LENGTHENED,
    NONE;

    public static JitterMatch checkJitter(
            final VariantReadContext readContext, final ReadContextMatcher matcher, final SAMRecord record, int readVarIndex)
    {
        if(readContext.AllRepeats.isEmpty())
            return JitterMatch.NONE;

        final byte[] readBases = record.getReadBases();
        final byte[] readQuals = record.getBaseQualities();
        boolean checkDoubleJitter = readContext.MaxRepeat != null && readContext.MaxRepeat.Count >= DOUBLE_JITTER_REPEAT_COUNT;

        // try each repeat covering the read context in turn
        for(RepeatInfo repeat : readContext.AllRepeats)
        {
            // test each jitter type in turn
            for(int i = 0; i <= 1; ++i)
            {
                JitterMatch jitterType = (i == 0) ? SHORTENED : LENGTHENED;

                for(int j = 0; j <= 1; ++j)
                {
                    // with or without awareness of an indel
                    boolean isIndelOffset = (j == 0);

                    for(int k = 0; k <= 1; ++k)
                    {
                        // before or after the variant
                        boolean jitterAtStart = (k == 0);

                        if(hasJitterMatchType(
                                repeat, readContext, matcher.altIndexLower(), matcher.altIndexUpper(), readVarIndex, readBases, readQuals,
                                jitterType, isIndelOffset, jitterAtStart, 1))
                        {
                            return jitterType;
                        }

                        if(isIndelOffset && checkDoubleJitter)
                        {
                            // for contexts with long repeats, test double jitter
                            if(hasJitterMatchType(
                                    repeat, readContext, matcher.altIndexLower(), matcher.altIndexUpper(), readVarIndex, readBases, readQuals,
                                    jitterType, isIndelOffset, jitterAtStart, 2))
                            {
                                return jitterType;
                            }
                        }
                    }
                }
            }
        }

        return JitterMatch.NONE;
    }

    @VisibleForTesting
    public static boolean hasJitterMatchType(
            final RepeatInfo repeat, final VariantReadContext readContext, int altIndexLower, int altIndexUpper,
            int readVarIndex, final byte[] readBases, final byte[] readQuals,
            final JitterMatch jitterType, boolean isIndelOffset, boolean jitterAtStart, int jitterCount)
    {
        int repeatLength = repeat.repeatLength() * jitterCount;
        int repeatEndIndex = repeat.endIndex();

        int readVarIndexOffset = readVarIndex - readContext.VarIndex;;
        int relativeRepeatIndex = repeat.Index + readVarIndexOffset;
        int relativeRepeatIndexEnd = repeatEndIndex + readVarIndexOffset;

        // cannot test a repeat if the read does not fully cover it
        if(relativeRepeatIndex < 0 || relativeRepeatIndexEnd >= readBases.length)
            return false;

        // for each read try shortening and then lengthening it
        int flankReadIndexStart = readVarIndex - readContext.leftLength();
        int flankReadIndexEnd = readVarIndex + readContext.rightLength() - 1;

        if(isIndelOffset && jitterAtStart)
        {
            // factor in repeats explained by an indel and which finish before the variant read index by shifting the implied flank start
            if(jitterType == SHORTENED)
                flankReadIndexStart += repeatLength;
            else
                flankReadIndexStart -= repeatLength;
        }

        // at least one flank must be fully present and the core must be covered
        int coreReadIndexStart = flankReadIndexStart + readContext.leftFlankLength();
        int coreReadIndexEnd = flankReadIndexEnd - readContext.rightFlankLength();

        if(coreReadIndexStart < 0 || coreReadIndexEnd >= readBases.length)
            return false;

        int readContextIndex = 0;
        int readIndex = flankReadIndexStart;
        boolean allMatched = true;
        boolean indexAdjusted = false;
        boolean altAdjusted = false;

        int altMatchCount = 0;
        int permittedLowQualMismatches = 1;
        int permittedLowQualRangeLower = readContext.VarIndex;
        int permittedLowQualRangeUpper = altIndexUpper + (jitterType == LENGTHENED ? 1 : 0);

        for(; readIndex <= flankReadIndexEnd; ++readIndex, ++readContextIndex)
        {
            if(readIndex < 0)
                continue;

            if(isIndelOffset && !indexAdjusted)
            {
                if(jitterType == SHORTENED && readContextIndex == repeatEndIndex - repeatLength)
                {
                    indexAdjusted = true;
                    readIndex -= repeatLength;
                }
                else if(jitterType == LENGTHENED && readContextIndex == repeatEndIndex + 1)
                {
                    indexAdjusted = true;
                    readContextIndex -= repeatLength;
                }

                if(readIndex < 0 || readContextIndex < 0)
                    return false;
            }
            else if(!isIndelOffset && !altAdjusted)
            {
                if(jitterAtStart)
                {
                    if((jitterType == SHORTENED && readContextIndex == repeat.Index)
                    || (jitterType == LENGTHENED && readContextIndex == repeat.Index - repeatLength))
                    {
                        altAdjusted = true;
                        altMatchCount = repeatLength;
                    }
                }
                else
                {
                    if((jitterType == SHORTENED && readContextIndex == repeatEndIndex)
                    || (jitterType == LENGTHENED && readContextIndex == repeatEndIndex + 1))
                    {
                        altAdjusted = true;
                        altMatchCount = repeatLength;
                    }
                }
            }

            if(readIndex >= readBases.length || readContextIndex >= readContext.ReadBases.length)
                break;

            byte readContextBase = readContext.ReadBases[readContextIndex];

            if(altMatchCount > 0)
            {
                --altMatchCount;

                if(jitterType == SHORTENED)
                {
                    // cannot skip testing the variant's bases itself
                    if(!readContext.variant().isIndel() && positionWithin(readContextIndex, altIndexLower, altIndexUpper))
                        return false;

                    continue;
                }

                // set to next expected repeat base
                int repeatBaseIndex = repeatLength - altMatchCount - 1;
                readContextBase = (byte)repeat.Bases.charAt(repeatBaseIndex);
            }

            if(readBases[readIndex] == readContextBase)
                continue;

            boolean differencePermitted = false;

            if(readQuals[readIndex] < MATCHING_BASE_QUALITY)
            {
                if(readContextIndex >= readContext.CoreIndexStart && readContextIndex <= readContext.CoreIndexEnd)
                {
                    // allow one mismatch in the core, but only outside a specific range from the variant index
                    if(permittedLowQualMismatches > 0)
                    {
                        if(readContextIndex < permittedLowQualRangeLower || readContextIndex > permittedLowQualRangeUpper)
                        {
                            --permittedLowQualMismatches;
                            differencePermitted = true;
                        }
                    }
                }
                else
                {
                    // low-qual outside the core is permitted
                    differencePermitted = true;
                }
            }

            if(!differencePermitted)
            {
                allMatched = false;
                break;
            }
        }

        return allMatched;
    }
}
