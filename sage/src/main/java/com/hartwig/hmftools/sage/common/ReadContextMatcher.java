package com.hartwig.hmftools.sage.common;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.sage.SageConstants.CORE_LOW_QUAL_MISMATCH_FACTOR;
import static com.hartwig.hmftools.sage.SageConstants.FLANK_LOW_QUAL_MISMATCHES;
import static com.hartwig.hmftools.sage.SageConstants.MATCHING_BASE_QUALITY;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.CORE;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.FULL;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.NONE;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.PARTIAL_CORE;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;

public class ReadContextMatcher
{
    private final VariantReadContext mContext;
    private final int mMaxCoreLowQualMatches;
    private final boolean mAllowWildcardMatchInCore;

    private final LowQualExclusion mLowQualExclusionRead;
    private final LowQualExclusion mLowQualExclusionRef;

    public static final byte WILDCARD_BASE = (byte) '.';

    public ReadContextMatcher(final VariantReadContext variantReadContext)
    {
        this(variantReadContext, true);
    }

    private class LowQualExclusion
    {
        public final int IndexLower;
        public final int IndexUpper;
        public final boolean IsRange;

        public LowQualExclusion(final int indexLower, final int indexUpper, final boolean isRange)
        {
            IndexLower = indexLower;
            IndexUpper = indexUpper;
            IsRange = isRange;
        }

        public boolean coversIndex(int index)
        {
            if(IsRange)
                return index >= IndexLower && index <= IndexUpper;
            else
                return index == IndexLower || index == IndexUpper;
        }

        public String toString() { return format("%d-%d %s", IndexLower, IndexUpper, IsRange ? "range" : "values"); }
    }

    public ReadContextMatcher(final VariantReadContext variantReadContext, boolean allowMismatches)
    {
        mContext = variantReadContext;

        mMaxCoreLowQualMatches = allowMismatches ? calcMaxLowQualCoreMismatches() : 0;
        mAllowWildcardMatchInCore = allowMismatches ? mContext.variant().isSNV() && !mContext.hasHomology() : false;

        if(mContext.variant().isIndel())
        {
            int altIndexLower = determineIndelLowQualLowerIndex(variantReadContext);
            int altIndexUpper = determineIndelLowQualUpperIndex(variantReadContext);

            mLowQualExclusionRead = new LowQualExclusion(altIndexLower, altIndexUpper, false);
            int refIndexDiff = mContext.leftFlankLength();
            mLowQualExclusionRef = new LowQualExclusion(altIndexLower - refIndexDiff, altIndexUpper - refIndexDiff, false);
        }
        else
        {
            // just the alt bases themselves - for both ref and read
            int altRange = mContext.variant().altLength() - 1;
            mLowQualExclusionRead = new LowQualExclusion(mContext.VarIndex, mContext.VarIndex + altRange, true);

            int refIndex = mContext.leftCoreLength();
            mLowQualExclusionRef = new LowQualExclusion(refIndex, refIndex + altRange, true);
        }
    }

    private static int determineIndelLowQualLowerIndex(final VariantReadContext readContext)
    {
        if(!readContext.variant().isInsert())
            return readContext.VarIndex;

        // return the last inserted base for insert, since this will be the first point of difference on the lower side
        return readContext.VarIndex + readContext.variant().indelLength();
    }

    private static int determineIndelLowQualUpperIndex(final VariantReadContext readContext)
    {
        // find the first base of difference (ref vs alt) up from the variant's position, and cap at the core indices
        final SimpleVariant variant = readContext.variant();

        int refIndex = readContext.variantRefIndex();
        int readIndex = readContext.VarIndex;

        for(; readIndex <= readContext.CoreIndexEnd & refIndex < readContext.RefBases.length; ++readIndex, ++refIndex)
        {
            if(readContext.RefBases[refIndex] != readContext.ReadBases[readIndex])
                break;
        }

        int upperRefIndex = variant.isInsert() ? readContext.VarIndex + variant.altLength() : readContext.VarIndex + 1;

        return max(min(readIndex, readContext.CoreIndexEnd), upperRefIndex);
    }

    private int calcMaxLowQualCoreMismatches()
    {
        int coreLength = mContext.coreLength();

        if(mContext.MaxRepeat != null)
        {
            // determine how many times the repeat bases are in the core
            int coreRepeatCount = 0;
            String coreStr = mContext.coreStr();
            int repeatIndex = coreStr.indexOf(mContext.MaxRepeat.Bases);

            while(repeatIndex >= 0)
            {
                ++coreRepeatCount;
                repeatIndex = coreStr.indexOf(mContext.MaxRepeat.Bases, repeatIndex + 1);
            }

            if(coreRepeatCount > 2)
            {
                int trimmedRepeat = (coreRepeatCount - 2) * mContext.MaxRepeat.repeatLength();
                coreLength = coreLength - trimmedRepeat;
            }
        }

        return calcMaxLowQualMismatches(coreLength);
    }

    private static int calcMaxLowQualMismatches(int segmentLength)
    {
        return 1 + (segmentLength / CORE_LOW_QUAL_MISMATCH_FACTOR);
    }

    public boolean coversVariant(final SAMRecord record, final int readVarIndex)
    {
        if(readVarIndex < 0)
            return false;

        int requiredReadIndexLower = readVarIndex + mContext.VarIndex - mContext.AltIndexLower;
        int requiredReadIndexUpper = readVarIndex + mContext.AltIndexUpper - mContext.VarIndex;

        // must cover from the first unambiguous ref vs alt bases on one side and the core in the opposite direction
        return requiredReadIndexLower >= 0 && requiredReadIndexUpper < record.getReadBases().length;
    }

    public ReadContextMatch determineReadMatch(final SAMRecord record, final int readVarIndex)
    {
        return determineReadMatch(record.getReadBases(), record.getBaseQualities(), readVarIndex, false);
    }

    public ReadContextMatch determineReadMatch(final byte[] readBases, final byte[] readQuals, final int readVarIndex, boolean skipRefMatch)
    {
        if(!skipRefMatch && coreMatchesRef(readBases, readQuals, readVarIndex))
            return ReadContextMatch.REF;

        ReadContextMatch coreMatch = determineCoreMatch(readBases, readQuals, readVarIndex);

        if(coreMatch == NONE)
            return NONE;

        BaseMatchType leftMatch = determineFlankMatch(readBases, readQuals, readVarIndex, true);
        BaseMatchType rightMatch = determineFlankMatch(readBases, readQuals, readVarIndex, false);

        if(rightMatch == BaseMatchType.MISMATCH || leftMatch == BaseMatchType.MISMATCH)
            return CORE;

        // flanks either matched or were incomplete
        if(coreMatch == PARTIAL_CORE)
        {
            if(leftMatch == BaseMatchType.MATCH || rightMatch == BaseMatchType.MATCH)
                return PARTIAL_CORE;
            else
                return CORE;
        }
        else
        {
            // again only one flank needs to be complete to be full (since no mismatches were found)
            if(leftMatch == BaseMatchType.MATCH || rightMatch == BaseMatchType.MATCH)
                return FULL;
            else
                return CORE;
        }
    }

    private boolean coreMatchesRef(final byte[] readBases, final byte[] readQuals, final int readVarIndex)
    {
        int refIndexStart = 0;

        int refIndex = mContext.variantRefIndex();
        int readIndexStart = readVarIndex - refIndex;
        int readIndexEnd = readIndexStart + mContext.RefBases.length - 1;

        int compareLength = mContext.RefBases.length;

        if(readIndexStart < 0 || readIndexEnd >= readBases.length)
            return false;

        return matches(
                mContext.RefBases, readBases, readQuals, refIndexStart, readIndexStart, compareLength,
                mMaxCoreLowQualMatches, mAllowWildcardMatchInCore, mLowQualExclusionRef);
    }

    private ReadContextMatch determineCoreMatch(final byte[] readBases, final byte[] readQuals, final int readVarIndex)
    {
        int readIndexStart = readVarIndex - mContext.leftCoreLength();
        int readIndexEnd = readVarIndex + mContext.rightCoreLength();

        // have already checked that the read covers the min alt range
        if(readIndexStart >= 0 && readIndexEnd < readBases.length)
        {
            if(matches(
                    mContext.ReadBases, readBases, readQuals, mContext.CoreIndexStart, readIndexStart, mContext.coreLength(),
                    mMaxCoreLowQualMatches, mAllowWildcardMatchInCore, mLowQualExclusionRead))
            {
                return CORE;
            }
            else
            {
                return NONE;
            }
        }
        else
        {
            int leftTrim = readIndexStart < 0 ? abs(readIndexStart) : 0;
            int rightTrim = readIndexEnd >= readBases.length ? abs(readBases.length - readIndexEnd - 1) : 0;

            readIndexStart += leftTrim;
            int coreIndexStart = mContext.CoreIndexStart + leftTrim;
            int compareLength = mContext.coreLength() - leftTrim - rightTrim;

            int maxLowQualMismatches = min(mMaxCoreLowQualMatches, calcMaxLowQualMismatches(compareLength));

            if(matches(
                    mContext.ReadBases, readBases, readQuals, coreIndexStart, readIndexStart, compareLength,
                    maxLowQualMismatches, mAllowWildcardMatchInCore, mLowQualExclusionRead))
            {
                return ReadContextMatch.PARTIAL_CORE;
            }
            else
            {
                return NONE;
            }
        }
    }

    private enum BaseMatchType
    {
        MATCH,
        INCOMPLETE,
        MISMATCH;
    }

    private BaseMatchType determineFlankMatch(final byte[] readBases, final byte[] readQuals, final int readVarIndex, boolean isLeft)
    {
        int readIndexStart, readIndexEnd, flankStartIndex, flankLength;

        if(isLeft)
        {
            readIndexStart = readVarIndex - mContext.leftLength();
            readIndexEnd = readVarIndex - mContext.leftCoreLength() - 1;
            flankStartIndex = 0;
            flankLength = mContext.leftFlankLength();
        }
        else
        {
            readIndexStart = readVarIndex + mContext.rightCoreLength() + 1;
            readIndexEnd = readVarIndex + mContext.rightLength() - 1;
            flankStartIndex = mContext.CoreIndexEnd + 1;
            flankLength = mContext.rightFlankLength();
        }

        boolean isPartial = false;

        if(readIndexStart < 0 || readIndexEnd >= readBases.length)
        {
            isPartial = true;

            int leftTrim = readIndexStart < 0 ? abs(readIndexStart) : 0;
            int rightTrim = readIndexEnd >= readBases.length ? abs(readBases.length - readIndexEnd - 1) : 0;

            flankLength -= leftTrim + rightTrim;

            readIndexStart += leftTrim;
            flankStartIndex += leftTrim;
        }

        if(flankLength <= 0)
            return BaseMatchType.INCOMPLETE;

        boolean flankMatch = matches(
                mContext.ReadBases, readBases, readQuals, flankStartIndex, readIndexStart, flankLength,
                FLANK_LOW_QUAL_MISMATCHES, false, mLowQualExclusionRead);

        if(!flankMatch)
            return BaseMatchType.MISMATCH;

        return isPartial ? BaseMatchType.INCOMPLETE : BaseMatchType.MATCH;
    }

    private static boolean matches(
            final byte[] bases, final byte[] readBases, final byte[] readQuals, final int baseIndexStart, final int readIndexStart,
            final int compareLength, final int maxLowQualMismatches, final boolean allowWildcardMismatches,
            final LowQualExclusion lowQualExclusion)
    {
        return matchType(
                bases, readBases, readQuals, baseIndexStart, readIndexStart, compareLength,
                maxLowQualMismatches, allowWildcardMismatches, lowQualExclusion) == BaseMatchType.MATCH;
    }

    private static BaseMatchType matchType(
            final byte[] bases, final byte[] readBases, @Nullable final byte[] readQuals, final int baseIndexStart, final int readIndexStart,
            final int compareLength, final int maxLowQualMismatches, final boolean allowWildcardMismatches,
            @Nullable final LowQualExclusion lowQualExclusion)
    {
        if(compareLength <= 0)
            return BaseMatchType.MISMATCH;

        // fail if asking to check bases beyond either array
        int baseIndexEnd = baseIndexStart + compareLength - 1;
        int readIndexEnd = readIndexStart + compareLength - 1;

        if(baseIndexEnd >= bases.length || readIndexEnd >= readBases.length)
            return BaseMatchType.INCOMPLETE;

        int mismatchCount = 0;

        for(int i = baseIndexStart, j = readIndexStart; i <= baseIndexEnd && j <= readIndexEnd; ++i, ++j)
        {
            if(bases[i] == readBases[j])
                continue;

            if(allowWildcardMismatches && readBases[j] == WILDCARD_BASE)
                continue;

            if(readQuals == null)
                return BaseMatchType.MISMATCH;

            if(lowQualExclusion != null && lowQualExclusion.coversIndex(i))
                return BaseMatchType.MISMATCH;

            if(readQuals[j] < MATCHING_BASE_QUALITY)
            {
                ++mismatchCount;

                if(mismatchCount > maxLowQualMismatches)
                    return BaseMatchType.MISMATCH;
            }
            else
            {
                return BaseMatchType.MISMATCH;
            }
        }

        return BaseMatchType.MATCH;
    }

    public static ReadContextMatch compareReadContexts(final VariantReadContext first, final VariantReadContext second)
    {
        // compare the core and flanks for the 2 contexts, not allowing for mismatches
        ReadContextMatcher matcher = new ReadContextMatcher(first, false);
        return matcher.determineReadMatch(second.ReadBases, null, second.VarIndex, true);
    }

    public double averageCoreQuality(final SAMRecord record, final int readVarIndex)
    {
        int readIndexStart = max(readVarIndex - mContext.leftCoreLength(), 0);
        int readIndexEnd = min(readVarIndex + mContext.rightCoreLength(), record.getReadBases().length - 1);

        int baseLength = readIndexEnd - readIndexStart + 1;

        if(baseLength <= 0)
            return 0;

        double quality = 0;

        for(int i = readIndexStart; i <= readIndexEnd; i++)
        {
            quality += record.getBaseQualities()[i];
        }

        return (int)round(quality / baseLength);
    }
}
