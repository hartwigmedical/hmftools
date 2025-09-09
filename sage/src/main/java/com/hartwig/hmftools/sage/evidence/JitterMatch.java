package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.sage.SageConstants.DOUBLE_JITTER_REPEAT_COUNT;
import static com.hartwig.hmftools.sage.SageConstants.MATCHING_BASE_QUALITY;
import static com.hartwig.hmftools.sage.quality.QualityCalculator.isMediumBaseQual;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.redux.BaseQualAdjustment;
import com.hartwig.hmftools.sage.common.ReadContextMatcher;
import com.hartwig.hmftools.sage.common.RepeatInfo;
import com.hartwig.hmftools.sage.common.VariantReadContext;

import htsjdk.samtools.SAMRecord;

public enum JitterMatch
{
    SHORTENED,
    LENGTHENED,
    NONE;

    public static boolean hasValidBaseQuals(
            final VariantReadContext readContext, final ReadContextMatcher matcher, final SAMRecord record, int readVarIndex)
    {
        int altIndexUpper = matcher.altIndexUpper();

        // for each read try shortening and then lengthening it
        int flankReadIndexStart = readVarIndex - readContext.leftLength();
        int flankReadIndexEnd = readVarIndex + readContext.rightLength() - 1;

        int coreReadIndexStart = flankReadIndexStart + readContext.leftFlankLength();
        int coreReadIndexEnd = flankReadIndexEnd - readContext.rightFlankLength();

        if(coreReadIndexStart < 0 || coreReadIndexEnd >= record.getBaseQualities().length)
            return false;

        int readContextIndex = 0;
        int readIndex = flankReadIndexStart;

        int permittedLowQualRangeLower = readContext.VarIndex;
        int permittedLowQualRangeUpper = altIndexUpper + 1;

        final byte[] readQuals = record.getBaseQualities();

        for(; readIndex <= flankReadIndexEnd; ++readIndex, ++readContextIndex)
        {
            boolean withinCore = readContextIndex >= readContext.CoreIndexStart && readContextIndex <= readContext.CoreIndexEnd;

            boolean withinPermittedRange = withinCore
                    && readContextIndex >= permittedLowQualRangeLower && readContextIndex <= permittedLowQualRangeUpper;

            if(withinPermittedRange && BaseQualAdjustment.isUncertainBaseFromQual(readQuals[readIndex]))
                return false;

            if(withinCore && isMediumBaseQual(readQuals[readIndex]))
                return false;
        }

        return true;
    }

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
                    // before or after the variant
                    boolean jitterAtStart = (j == 0);

                    if(hasJitterMatchType(
                            repeat, readContext, matcher.altIndexLower(), matcher.altIndexUpper(), readVarIndex, readBases, readQuals,
                            jitterType, jitterAtStart, 1))
                    {
                        return jitterType;
                    }

                    if(checkDoubleJitter)
                    {
                        // for contexts with long repeats, test double jitter
                        if(hasJitterMatchType(
                                repeat, readContext, matcher.altIndexLower(), matcher.altIndexUpper(), readVarIndex, readBases, readQuals,
                                jitterType, jitterAtStart, 2))
                        {
                            return jitterType;
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
            final JitterMatch jitterType, boolean jitterAtStart, int jitterCount)
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

        if(jitterAtStart)
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

        int altMatchCount = 0;
        int permittedLowQualMismatches = 1;
        int permittedLowQualRangeLower = readContext.VarIndex;
        int permittedLowQualRangeUpper = altIndexUpper + (jitterType == LENGTHENED ? 1 : 0);

        for(; readIndex <= flankReadIndexEnd; ++readIndex, ++readContextIndex)
        {
            if(readIndex < 0)
                continue;

            if(!indexAdjusted)
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

            boolean withinCore = readContextIndex >= readContext.CoreIndexStart && readContextIndex <= readContext.CoreIndexEnd;

            boolean withinPermittedRange = withinCore
                    && readContextIndex >= permittedLowQualRangeLower && readContextIndex <= permittedLowQualRangeUpper;

            if(withinPermittedRange && BaseQualAdjustment.isUncertainBaseFromQual(readQuals[readIndex]))
                return false;

            if(withinCore && isMediumBaseQual(readQuals[readIndex]))
                return false;

            if(readQuals[readIndex] < MATCHING_BASE_QUALITY)
            {
                if(withinCore)
                {
                    // allow one mismatch in the core, but only outside a specific range from the variant index
                    if(permittedLowQualMismatches > 0)
                    {
                        if(!withinPermittedRange)
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
