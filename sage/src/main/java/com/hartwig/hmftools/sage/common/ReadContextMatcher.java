package com.hartwig.hmftools.sage.common;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.sage.SageConstants.CORE_LOW_QUAL_MISMATCH_FACTOR;
import static com.hartwig.hmftools.sage.SageConstants.FLANK_LOW_QUAL_MISMATCHES;
import static com.hartwig.hmftools.sage.SageConstants.LONG_GERMLINE_INSERT_READ_VS_REF_DIFF;
import static com.hartwig.hmftools.sage.SageConstants.MATCHING_BASE_QUALITY;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.CORE;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.FULL;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.NONE;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.PARTIAL_CORE;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.SIMPLE_ALT;
import static com.hartwig.hmftools.sage.common.SimpleVariant.isLongInsert;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class ReadContextMatcher
{
    private final VariantReadContext mContext;
    private final int mMaxCoreLowQualMatches;

    private final boolean mIsReference;
    private final int mAltIndexLower; // the first lower base relative to the index where the ref and alt differ
    private final int mAltIndexUpper;

    private final LowQualExclusion mLowQualExclusionRead;
    private final LowQualExclusion mLowQualExclusionRef;

    private class LowQualExclusion
    {
        public final List<Integer> Indices;
        public final boolean IsRange;

        private int mIndexOffset;

        public LowQualExclusion(final int indexLower, final int indexUpper)
        {
            Indices = List.of(indexLower, indexUpper);
            IsRange = true;
            mIndexOffset = 0;
        }

        public LowQualExclusion(final List<Integer> indices)
        {
            Indices = indices;
            IsRange = false;
            mIndexOffset = 0;
        }

        public void setIndexOffset(int offset) { mIndexOffset = offset; }

        public boolean coversIndex(int index)
        {
            index += mIndexOffset;

            if(IsRange)
                return index >= Indices.get(0) && index <= Indices.get(1);
            else
                return Indices.contains(index);
        }

        public String toString() { return format("%s: %s", IsRange ? "range" : "values", Indices); }
    }

    public ReadContextMatcher(final VariantReadContext readContext, boolean allowMismatches, boolean isReference)
    {
        mContext = readContext;

        mMaxCoreLowQualMatches = allowMismatches ? calcMaxLowQualCoreMismatches() : 0;

        mIsReference = isReference;

        int altIndexLower = readContext.VarIndex;
        int altIndexUpper = determineAltIndexUpper(readContext.variant(), readContext.VarIndex, readContext.Homology);

        if(!mIsReference && !readContext.AllRepeats.isEmpty())
        {
            // expand alt boundaries to cover any repeats on that side
            int minRepeatIndex = readContext.AllRepeats.get(0).Index - 1;
            int maxRepeatIndex = readContext.AllRepeats.get(0).endIndex() + 1;

            for(int i = 1; i < readContext.AllRepeats.size(); ++i)
            {
                minRepeatIndex = min(minRepeatIndex, readContext.AllRepeats.get(i).Index - 1);
                maxRepeatIndex = max(maxRepeatIndex, readContext.AllRepeats.get(i).endIndex() + 1);
            }

            altIndexLower = min(altIndexLower, minRepeatIndex);
            altIndexUpper = max(altIndexUpper, maxRepeatIndex);
        }

        mAltIndexLower = altIndexLower;
        mAltIndexUpper = altIndexUpper;

        if(allowMismatches)
        {
            if(mContext.variant().isIndel())
            {
                Set<Integer> excludedBases = Sets.newHashSet();
                int lowQualIndexLower = determineIndelLowQualLowerIndex(readContext);
                int lowQualIndexUpper = determineIndelLowQualRefReadDiffIndex(readContext);
                excludedBases.add(lowQualIndexLower);
                excludedBases.add(lowQualIndexUpper);

                if(mContext.Homology != null)
                {
                    excludedBases.add(mContext.VarIndex);
                    excludedBases.add(mContext.VarIndex + 1);
                    excludedBases.add(mAltIndexUpper - 1);
                    excludedBases.add(mAltIndexUpper);
                }

                mLowQualExclusionRead = new LowQualExclusion(excludedBases.stream().collect(Collectors.toList()));

                int refIndexDiff = mContext.leftFlankLength();
                mLowQualExclusionRef = new LowQualExclusion(List.of(lowQualIndexLower - refIndexDiff, lowQualIndexUpper - refIndexDiff));
            }
            else
            {
                // just the alt bases themselves - for both ref and read
                int altRange = mContext.variant().altLength() - 1;
                mLowQualExclusionRead = new LowQualExclusion(mContext.VarIndex, mContext.VarIndex + altRange);

                int refIndex = mContext.leftCoreLength();
                mLowQualExclusionRef = new LowQualExclusion(refIndex, refIndex + altRange);
            }
        }
        else
        {
            mLowQualExclusionRead = null;
            mLowQualExclusionRef = null;
        }
    }

    public boolean isReference() { return mIsReference; }
    public int altIndexLower() { return mAltIndexLower; }
    public int altIndexUpper() { return mAltIndexUpper; }

    public void setRealignmentIndexOffset(int offset) { mLowQualExclusionRead.setIndexOffset(offset); }
    public void clearRealignmentIndexOffset() { mLowQualExclusionRead.setIndexOffset(0); }

    private static int determineAltIndexUpper(final SimpleVariant variant, final int readVarIndex, final Microhomology homology)
    {
        if(variant.isInsert())
            return readVarIndex + (homology != null ? homology.Length : 0) + 1 + abs(variant.indelLength());
        else if(variant.isDelete())
            return readVarIndex + (homology != null ? homology.Length : 0) + 1;
        else
            return readVarIndex + variant.altLength() - 1;
    }

    private static int determineIndelLowQualLowerIndex(final VariantReadContext readContext)
    {
        if(!readContext.variant().isInsert())
            return readContext.VarIndex;

        // return the last inserted base for insert, since this will be the first point of difference on the lower side
        return readContext.VarIndex + readContext.variant().indelLength();
    }

    private static int determineIndelLowQualRefReadDiffIndex(final VariantReadContext readContext)
    {
        // find the first base of difference (ref vs alt) up from the variant's position, and cap at the core indices
        int refIndex = readContext.variantRefIndex();
        int readIndex = readContext.VarIndex;

        for(; readIndex <= readContext.CoreIndexEnd & refIndex < readContext.RefBases.length; ++readIndex, ++refIndex)
        {
            if(readContext.RefBases[refIndex] != readContext.ReadBases[readIndex])
                break;
        }

        return min(readIndex, readContext.CoreIndexEnd);
    }

    private int calcMaxLowQualCoreMismatches()
    {
        int coreLength = mContext.coreLength();

        if(mContext.MaxRepeat != null)
        {
            // determine how many times the repeat bases are in the core, starting from the repeat itself
            int coreRepeatCount = 0;
            String coreStr = mContext.coreStr();
            int repeatIndexStart = mContext.MaxRepeat.Index - mContext.leftFlankLength();
            int repeatIndex = coreStr.indexOf(mContext.MaxRepeat.Bases, repeatIndexStart);

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

    public boolean coversVariant(final SAMRecord record, final int readVarIndex) { return coversVariant(record.getReadBases(), readVarIndex); }

    public boolean coversVariant(final byte[] readBases, final int readVarIndex)
    {
        if(readVarIndex < 0)
            return false;

        int requiredReadIndexLower = readVarIndex - mContext.VarIndex + mAltIndexLower;
        int requiredReadIndexUpper = readVarIndex - mContext.VarIndex + mAltIndexUpper;

        // must cover from the first unambiguous ref vs alt bases on one side and the core in the opposite direction
        return requiredReadIndexLower >= 0 && requiredReadIndexUpper < readBases.length;
    }

    public ReadContextMatch determineReadMatch(final SAMRecord record, final int readVarIndex)
    {
        ReadContextMatch matchType = determineReadMatch(record.getReadBases(), record.getBaseQualities(), readVarIndex, false);

        if(matchType == NONE && mIsReference && checkSimpleAltMatch(record, readVarIndex))
            return SIMPLE_ALT;

        return matchType;
    }

    public ReadContextMatch determineReadMatch(final byte[] readBases, final byte[] readQuals, final int readVarIndex, boolean skipRefMatch)
    {
        if(!skipRefMatch && coreMatchesRef(readBases, readQuals, readVarIndex))
            return ReadContextMatch.REF;

        ReadContextMatch coreMatch = determineCoreMatch(readBases, readQuals, readVarIndex);

        if(coreMatch == NONE)
        {
            if(mIsReference && isLongInsert(mContext.variant()))
            {
                String extendedRefBases = mContext.refBases() + mContext.extendedRefBases();

                int[] refReadMatches = countRefAndReadDifferences(
                        readBases, mContext.ReadBases, extendedRefBases.getBytes(), readVarIndex, mContext.VarIndex, mContext.variantRefIndex());

                if(refReadMatches[1] - refReadMatches[0] >= LONG_GERMLINE_INSERT_READ_VS_REF_DIFF)
                    return PARTIAL_CORE;
            }

            return NONE;
        }

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

        // allow partial ref core to be compared since, where required, checks have already been made that the core is sufficiently covered
        if(readIndexStart < 0)
        {
            refIndexStart += abs(readIndexStart);
            readIndexStart = 0;
        }

        readIndexEnd = min(readIndexEnd, readBases.length - 1);

        int compareLength = readIndexEnd - readIndexStart + 1;

        return matches(
                mContext.RefBases, readBases, readQuals, refIndexStart, readIndexStart, compareLength,
                mMaxCoreLowQualMatches, mLowQualExclusionRef);
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
                    mMaxCoreLowQualMatches, mLowQualExclusionRead))
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
                    maxLowQualMismatches, mLowQualExclusionRead))
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
                FLANK_LOW_QUAL_MISMATCHES, mLowQualExclusionRead);

        if(!flankMatch)
            return BaseMatchType.MISMATCH;

        return isPartial ? BaseMatchType.INCOMPLETE : BaseMatchType.MATCH;
    }

    private static boolean matches(
            final byte[] bases, final byte[] readBases, final byte[] readQuals, final int baseIndexStart, final int readIndexStart,
            final int compareLength, final int maxLowQualMismatches, final LowQualExclusion lowQualExclusion)
    {
        return matchType(
                bases, readBases, readQuals, baseIndexStart, readIndexStart, compareLength,
                maxLowQualMismatches, lowQualExclusion) == BaseMatchType.MATCH;
    }

    private static BaseMatchType matchType(
            final byte[] bases, final byte[] readBases, @Nullable final byte[] readQuals, final int baseIndexStart, final int readIndexStart,
            final int compareLength, final int maxLowQualMismatches, @Nullable final LowQualExclusion lowQualExclusion)
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

    private boolean checkSimpleAltMatch(final SAMRecord record, final int readVarIndex)
    {
        SimpleVariant variant = mContext.variant();

        if(variant.isIndel())
        {
            // check that the indel exists at this location and matches
            int readIndex = 0;
            for(CigarElement element : record.getCigar().getCigarElements())
            {
                if(readIndex == readVarIndex + 1)
                {
                    if(element.getLength() != variant.indelLengthAbs())
                        return false;

                    if(variant.isDelete())
                    {
                        return element.getOperator() == D;
                    }
                    else
                    {
                        if(element.getOperator() != I)
                            return false;

                        // check that inserted bases match
                        for(int i = 0; i < element.getLength(); ++i)
                        {
                            byte altBase = (byte)variant.alt().charAt(i + 1);
                            byte readBase = record.getReadBases()[readIndex + i];

                            if(readBase != altBase)
                                return false;
                        }

                        return true;
                    }
                }

                if(element.getOperator().consumesReadBases())
                    readIndex += element.getLength();
            }

            return false;
        }
        else
        {
            for(int i = 0 ; i < variant.altLength(); ++i)
            {
                byte altBase = (byte)variant.alt().charAt(i);
                byte readBase = record.getReadBases()[readVarIndex + i];

                if(readBase != altBase)
                    return false;
            }

            return true;
        }
    }

    private static int[] countRefAndReadDifferences(
            final byte[] readBases, final byte[] readContextBases, final byte[] refBases,
            final int readIndexStart, final int readContextIndexStart, final int refIndexStart)
    {
        // fail if asking to check bases beyond either array
        int readContextMatches = 0;
        int refMatches = 0;

        int i = 0;

        while(true)
        {
            int readBaseIndex = readIndexStart + i;

            if(readBaseIndex >= readBases.length)
                break;

            int readContextIndex = readContextIndexStart + i;

            if(readContextIndex >= readContextBases.length)
                break;

            int refIndex = refIndexStart + i;

            if(refIndex >= refBases.length)
                break;

            if(readBases[readBaseIndex] == readContextBases[readContextIndex])
                ++readContextMatches;

            if(readBases[readBaseIndex] == refBases[refIndex])
                ++refMatches;

            ++i;
        }

        return new int[] { refMatches, readContextMatches };
    }
}
