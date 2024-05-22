package com.hartwig.hmftools.sage.common;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

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

    private final int[] mLowQualExclusionRangeRead; // for SNV & MNVs, the indices covering the variant & excluding low-qual mismatches
    private final int[] mLowQualExclusionRangeRef;

    public static final byte WILDCARD_BASE = (byte) '.';

    public ReadContextMatcher(final VariantReadContext variantReadContext)
    {
        this(variantReadContext, true);
    }

    public ReadContextMatcher(final VariantReadContext variantReadContext, boolean allowMismatches)
    {
        mContext = variantReadContext;

        mMaxCoreLowQualMatches = allowMismatches ? calcMaxLowQualCoreMismatches() : 0;
        mAllowWildcardMatchInCore = allowMismatches ? mContext.variant().isSNV() && !mContext.hasHomology() : false;

        if(mContext.variant().isIndel())
        {
            mLowQualExclusionRangeRead = new int[] { mContext.AltIndexUpper, mContext.AltIndexUpper };
            mLowQualExclusionRangeRef = mLowQualExclusionRangeRead;
        }
        else
        {
            // just the alt bases themselves - for both ref and read
            int altLength = mContext.variant().altLength();
            mLowQualExclusionRangeRead = new int[] { mContext.VarReadIndex, mContext.VarReadIndex + altLength - 1 };
            int refIndex = mContext.leftCoreLength();
            mLowQualExclusionRangeRef = new int[] { refIndex, refIndex + altLength - 1 };
        }
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

    public boolean coversVariant(final SAMRecord record, final int readIndex)
    {
        int requiredReadIndexLower = readIndex + mContext.VarReadIndex - mContext.AltIndexLower;
        int requiredReadIndexUpper = readIndex + mContext.AltIndexUpper - mContext.VarReadIndex;

        // must cover from the first unambiguous ref vs alt bases on one side and the core in the opposite direction
        return requiredReadIndexLower >= 0 && requiredReadIndexUpper < record.getReadBases().length;
    }

    public ReadContextMatch determineReadMatch(final SAMRecord record, final int readIndex)
    {
        return determineReadMatch(record.getReadBases(), record.getBaseQualities(), readIndex, false);
    }

    public ReadContextMatch determineReadMatch(final byte[] readBases, final byte[] readQuals, final int readIndex, boolean skipRefMatch)
    {
        if(!skipRefMatch && coreMatchesRef(readBases, readQuals, readIndex))
            return ReadContextMatch.REF;

        ReadContextMatch coreMatch = determineCoreMatch(readBases, readQuals, readIndex);

        if(coreMatch == NONE)
            return NONE;

        BaseMatchType leftMatch = determineFlankMatch(readBases, readQuals, readIndex, true);
        BaseMatchType rightMatch = determineFlankMatch(readBases, readQuals, readIndex, false);

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

    private boolean coreMatchesRef(final byte[] readBases, final byte[] readQuals, final int readIndex)
    {
        int refIndexStart = 0;

        int readIndexStart = readIndex - mContext.leftCoreLength();
        int readIndexEnd = readIndexStart + mContext.RefBases.length - 1;

        int compareLength = mContext.RefBases.length;

        if(readIndexStart < 0 || readIndexEnd >= readBases.length)
            return false;

        return matches(
                mContext.RefBases, readBases, readQuals, refIndexStart, readIndexStart, compareLength,
                mMaxCoreLowQualMatches, mAllowWildcardMatchInCore, mLowQualExclusionRangeRef);
    }

    private ReadContextMatch determineCoreMatch(final byte[] readBases, final byte[] readQuals, final int readIndex)
    {
        int readIndexStart = readIndex - mContext.leftCoreLength();
        int readIndexEnd = readIndex + mContext.rightCoreLength();

        // have already checked that the read covers the min alt range
        if(readIndexStart >= 0 && readIndexEnd < readBases.length)
        {
            if(matches(
                    mContext.ReadBases, readBases, readQuals, mContext.CoreIndexStart, readIndexStart, mContext.coreLength(),
                    mMaxCoreLowQualMatches, mAllowWildcardMatchInCore, mLowQualExclusionRangeRead))
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
                    maxLowQualMismatches, mAllowWildcardMatchInCore, mLowQualExclusionRangeRead))
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

    private BaseMatchType determineFlankMatch(final byte[] readBases, final byte[] readQuals, final int readIndex, boolean isLeft)
    {
        int readIndexStart, readIndexEnd, flankStartIndex, flankLength;

        if(isLeft)
        {
            readIndexStart = readIndex - mContext.leftLength();
            readIndexEnd = readIndex - mContext.leftCoreLength() - 1;
            flankStartIndex = 0;
            flankLength = mContext.leftFlankLength();
        }
        else
        {
            readIndexStart = readIndex + mContext.rightCoreLength() + 1;
            readIndexEnd = readIndex + mContext.rightLength() - 1;
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
                FLANK_LOW_QUAL_MISMATCHES, false, mLowQualExclusionRangeRead);

        if(!flankMatch)
            return BaseMatchType.MISMATCH;

        return isPartial ? BaseMatchType.INCOMPLETE : BaseMatchType.MATCH;
    }

    private static boolean matches(
            final byte[] bases, final byte[] readBases, final byte[] readQuals, final int baseIndexStart, final int readIndexStart,
            final int compareLength, final int maxLowQualMismatches, final boolean allowWildcardMismatches,
            final int[] lowQualExclusionRange)
    {
        return matchType(
                bases, readBases, readQuals, baseIndexStart, readIndexStart, compareLength,
                maxLowQualMismatches, allowWildcardMismatches, lowQualExclusionRange) == BaseMatchType.MATCH;
    }

    private static BaseMatchType matchType(
            final byte[] bases, final byte[] readBases, @Nullable final byte[] readQuals, final int baseIndexStart, final int readIndexStart,
            final int compareLength, final int maxLowQualMismatches, final boolean allowWildcardMismatches,
            @Nullable final int[] lowQualExclusionRange)
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

            if(lowQualExclusionRange != null && i >= lowQualExclusionRange[0] && i <= lowQualExclusionRange[1])
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
        return matcher.determineReadMatch(second.ReadBases, null, second.VarReadIndex, true);
    }

    public double averageCoreQuality(final SAMRecord record, final int readIndex)
    {
        int readIndexStart = max(readIndex - mContext.leftCoreLength(), 0);
        int readIndexEnd = min(readIndex + mContext.rightCoreLength(), record.getReadBases().length - 1);

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
