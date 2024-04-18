package com.hartwig.hmftools.sage.common;

import static java.lang.Math.abs;
import static java.lang.Math.round;

import static com.hartwig.hmftools.sage.SageConstants.CORE_LOW_QUAL_MISMATCH_FACTOR;
import static com.hartwig.hmftools.sage.SageConstants.MATCHING_BASE_QUALITY;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.CORE;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.FULL;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.NONE;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.PARTIAL_CORE;

import com.google.common.annotations.VisibleForTesting;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;

public class ReadContextMatcher
{
    private final VariantReadContext mContext;
    private final int mMaxCoreLowQualMatches;
    private final boolean mAllowWildcardMatchInCore;

    private final int[] mLowQualExclusionRange;

    public static final byte WILDCARD_BASE = (byte) '.';

    public ReadContextMatcher(final VariantReadContext variantReadContext)
    {
        this(variantReadContext, true);
    }

    public ReadContextMatcher(final VariantReadContext variantReadContext, boolean allowMismatches)
    {
        mContext = variantReadContext;

        mMaxCoreLowQualMatches = allowMismatches ? calcMaxLowQualMismatches(mContext.coreLength()) : 0;
        mAllowWildcardMatchInCore = allowMismatches ? mContext.variant().isSNV() && !mContext.hasHomology() : false;

        if(mContext.variant().isIndel())
        {
            mLowQualExclusionRange = null;
        }
        else
        {
            mLowQualExclusionRange = new int[] { mContext.VarReadIndex, mContext.VarReadIndex + mContext.variant().alt().length() - 1 };
        }
    }

    private static int calcMaxLowQualMismatches(int segmentLength)
    {
        return 1 + (segmentLength / CORE_LOW_QUAL_MISMATCH_FACTOR);
    }

    public boolean coversVariant(final SAMRecord record, final int readIndex)
    {
        // must cover from the first unambiguous ref vs alt bases on one side and the core in the opposite direction
        int readLowerExtensionLength = readIndex;
        int readUpperExtensionLength = record.getReadBases().length - readIndex - 1;

        if(readLowerExtensionLength < mContext.VarReadIndex - mContext.AltIndexLower)
            return false;
        else if(readUpperExtensionLength < mContext.AltIndexUpper - mContext.VarReadIndex)
            return false;
        else
            return true;
    }

    public ReadContextMatch determineReadMatch(final SAMRecord record, final int readIndex)
    {
        return determineReadMatch(record.getReadBases(), record.getBaseQualities(), readIndex);
    }

    public ReadContextMatch determineReadMatch(final byte[] readBases, final byte[] readQuals, final int readIndex)
    {
        if(coreMatchesRef(readBases, readQuals, readIndex))
            return ReadContextMatch.REF;

        ReadContextMatch coreMatch = determineCoreMatch(readBases, readQuals, readIndex);

        if(coreMatch == NONE)
            return NONE;

        boolean leftMatch = determineFlankMatch(readBases, readQuals, readIndex, true);
        boolean rightMatch = determineFlankMatch(readBases, readQuals, readIndex, false);

        if(coreMatch == PARTIAL_CORE && (leftMatch || rightMatch))
            return PARTIAL_CORE;

        return (coreMatch == CORE && leftMatch && rightMatch) ? FULL : CORE;
    }

    private boolean coreMatchesRef(final byte[] readBases, final byte[] readQuals, final int readIndex)
    {
        int refIndex = mContext.refIndex();
        int refCoreIndexStart = refIndex - mContext.leftCoreLength();

        int readIndexStart = readIndex - mContext.leftCoreLength();
        int readIndexEnd = readIndex + mContext.rightCoreLength();

        if(readIndexStart < 0 || readIndexEnd >= readBases.length)
            return false;

        return matches(
                mContext.RefBases, readBases, readQuals, refCoreIndexStart, readIndexStart, mContext.coreLength(),
                mMaxCoreLowQualMatches, mAllowWildcardMatchInCore, mLowQualExclusionRange);
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
                    mMaxCoreLowQualMatches, mAllowWildcardMatchInCore, mLowQualExclusionRange))
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
            int rightTrim = readIndexEnd >= readBases.length ? abs(readBases.length - readIndexEnd + 1) : 0;

            readIndexStart += leftTrim;
            int coreIndexStart = mContext.CoreIndexStart + leftTrim;
            int compareLength = mContext.coreLength() - leftTrim - rightTrim;

            int maxLowQualMismatches = calcMaxLowQualMismatches(compareLength);

            if(matches(
                    mContext.ReadBases, readBases, readQuals, coreIndexStart, readIndexStart, compareLength,
                    maxLowQualMismatches, mAllowWildcardMatchInCore, mLowQualExclusionRange))
            {
                return ReadContextMatch.PARTIAL_CORE;
            }
            else
            {
                return NONE;
            }
        }
    }

    private static final int ANY_LOW_BASE_MISMATCH = -1;

    private boolean determineFlankMatch(final byte[] readBases, final byte[] readQuals, final int readIndex, boolean isLeft)
    {
        boolean isPartial = false;
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

        if(readIndexStart < 0 || readIndexEnd >= readBases.length)
        {
            int leftTrim = readIndexStart < 0 ? abs(readIndexStart) : 0;
            int rightTrim = readIndexEnd >= readBases.length ? abs(readBases.length - readIndexEnd + 1) : 0;

            flankLength -= leftTrim + rightTrim;
            isPartial = true;

            readIndexStart += leftTrim;
            flankStartIndex += leftTrim;
        }

        if(flankLength <= 0)
            return false;

        return matches(
                mContext.ReadBases, readBases, readQuals, flankStartIndex, readIndexStart, flankLength,
                ANY_LOW_BASE_MISMATCH, false, mLowQualExclusionRange);
    }

    public static boolean matches(
            final byte[] bases, final byte[] readBases, @Nullable final byte[] readQuals, final int baseIndexStart, final int readIndexStart,
            final int compareLength, final int maxLowQualMismatches, final boolean allowWildcardMismatches,
            @Nullable final int[] lowQualExclusionRange)
    {
        if(compareLength <= 0)
            return false;

        int mismatchCount = 0;

        for(int i = baseIndexStart, j = readIndexStart; i < baseIndexStart + compareLength && j < readIndexStart + compareLength; ++i, ++j)
        {
            if(bases[i] == readBases[j])
                continue;

            if(allowWildcardMismatches && readBases[j] == WILDCARD_BASE)
                continue;

            if(readQuals == null)
                return false;

            if(lowQualExclusionRange != null && i >= lowQualExclusionRange[0] && i <= lowQualExclusionRange[1])
                return false;

            if(readQuals[j] < MATCHING_BASE_QUALITY)
            {
                if(maxLowQualMismatches == ANY_LOW_BASE_MISMATCH)
                    continue;

                ++mismatchCount;

                if(mismatchCount > maxLowQualMismatches)
                    return false;
            }
            else
            {
                return false;
            }
        }

        return true;
    }

    public static ReadContextMatch compareReadContexts(final VariantReadContext first, final VariantReadContext second)
    {
        // compare the core and flanks for the 2 contexts, not allowing for mismatches
        ReadContextMatcher matcher = new ReadContextMatcher(first, false);
        return matcher.determineReadMatch(second.ReadBases, null, second.VarReadIndex);
    }

    @VisibleForTesting
    public ReadContextMatch compareReadContexts(final VariantReadContext other)
    {
        return determineReadMatch(other.ReadBases, null, other.VarReadIndex);
    }

    public double averageCoreQuality(final SAMRecord record, final int readIndex)
    {
        // CLEAN-UP: soon to be replaced with new MSI model??
        int readIndexStart = readIndex - mContext.leftCoreLength();
        int readIndexEnd = readIndex + mContext.rightCoreLength();

        if(readIndexStart < 0 || readIndexEnd >= record.getReadBases().length)
            return 0;

        double quality = 0;
        int baseLength = readIndexEnd - readIndexStart + 1;

        for(int i = readIndexStart; i <= readIndexEnd; i++)
        {
            quality += record.getBaseQualities()[i];
        }

        return (int)round(quality / baseLength);
    }
}
