package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;

import javax.annotation.Nullable;

import com.hartwig.hmftools.sage.common.VariantReadContext;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class ReadContextClassifier
{
    // TODO: Update this.
    public final static int HIGH_BASE_QUAL_CUTOFF = 37;

    private final VariantReadContext mVariantReadContext;

    public ReadContextClassifier(final VariantReadContext variantReadContext)
    {
        mVariantReadContext = variantReadContext;
    }

    // TODO: What happens if we go into soft-clips?
    // TODO: Other cigar types?
    // TODO: Fragment data.
    @Nullable
    public ReadContextCounter.MatchType classifyRead(final SAMRecord read)
    {
        if(read.getReadUnmappedFlag())
        {
            return null;
        }

        // TODO: Chromosome matches?

        if(!positionsOverlap(mVariantReadContext.AlignmentStart, mVariantReadContext.AlignmentEnd, read.getAlignmentStart(), read.getAlignmentEnd()))
        {
            return null;
        }

        CoreMatchType coreMatchType = coreMatch(read);
        if(coreMatchType == null)
        {
            return ReadContextCounter.MatchType.NONE;
        }

        boolean matchesLeftFlank = isLeftFlankMatch(read);
        boolean matchesRightFlank = isRightFlankMatch(read);
        ReadContextCounter.MatchType matchType = getMatchType(coreMatchType.Extent, matchesLeftFlank, matchesRightFlank);

        if(!coreMatchType.MatchesRef)
        {
            return matchType;
        }

        return matchType == ReadContextCounter.MatchType.NONE ? ReadContextCounter.MatchType.NONE : ReadContextCounter.MatchType.REF;
    }

    @NotNull
    private static ReadContextCounter.MatchType getMatchType(final CoreMatchExtent coreMatchExtent, final boolean matchesLeftFlank,
            final boolean matchesRightFlank)
    {
        if(coreMatchExtent == CoreMatchExtent.FULL && (matchesLeftFlank || matchesRightFlank))
        {
            return ReadContextCounter.MatchType.FULL;
        }

        if(matchesLeftFlank && coreMatchExtent == CoreMatchExtent.LEFT_PARTIAL)
        {
            return ReadContextCounter.MatchType.PARTIAL;
        }

        if(matchesRightFlank && coreMatchExtent == CoreMatchExtent.RIGHT_PARTIAL)
        {
            return ReadContextCounter.MatchType.PARTIAL;
        }

        if(coreMatchExtent == CoreMatchExtent.FULL)
        {
            return ReadContextCounter.MatchType.CORE;
        }

        return ReadContextCounter.MatchType.NONE;
    }

    // TODO: Do we need ALT_FULL?
    // TODO: Rename?
    private enum CoreMatchExtent
    {
        LEFT_PARTIAL,
        RIGHT_PARTIAL,
        FULL,
    }

    private static class CoreMatchType
    {
        public final CoreMatchExtent Extent;
        public final boolean MatchesRef;

        public CoreMatchType(final CoreMatchExtent extent, boolean matchesRef)
        {
            Extent = extent;
            MatchesRef = matchesRef;
        }
    }

    // TODO: Refactor and simplify.
    @Nullable
    private CoreMatchType coreMatch(final SAMRecord read)
    {
        String readString = read.getReadString();
        String coreString = mVariantReadContext.coreStr();

        // check variant
        String ref = mVariantReadContext.ref();
        String alt = mVariantReadContext.alt();
        int variantStartPos = mVariantReadContext.AlignmentStart + mVariantReadContext.VarReadIndex;
        String readAtVariant =
                readString.substring(variantStartPos - read.getUnclippedStart(), variantStartPos - read.getUnclippedStart() + alt.length());
        boolean matchesRef;
        if(readAtVariant.equals(ref))
        {
            matchesRef = true;
        }
        else if(readAtVariant.equals(alt))
        {
            matchesRef = false;
        }
        else
        {
            return null;
        }

        // check left of variant
        boolean matchesLeft = true;
        int coreStartPos = mVariantReadContext.AlignmentStart + mVariantReadContext.CoreIndexStart;
        for(int i = mVariantReadContext.CoreIndexStart; i < mVariantReadContext.VarReadIndex; ++i)
        {
            int coreIndex = i - mVariantReadContext.CoreIndexStart;
            int readIndex = coreIndex + coreStartPos - read.getUnclippedStart();
            if(coreString.charAt(coreIndex) != readString.charAt(readIndex))
            {
                matchesLeft = false;
                break;
            }
        }

        // check right of variant
        boolean matchesRight = true;
        for(int i = mVariantReadContext.VarReadIndex + alt.length(); i <= mVariantReadContext.CoreIndexEnd; ++i)
        {
            int coreIndex = i - mVariantReadContext.CoreIndexStart;
            int readIndex = coreIndex + coreStartPos - read.getUnclippedStart();
            if(coreString.charAt(coreIndex) != readString.charAt(readIndex))
            {
                matchesRight = false;
                break;
            }
        }

        if(matchesLeft && matchesRight)
        {
            return new CoreMatchType(CoreMatchExtent.FULL, matchesRef);
        }

        if(matchesLeft)
        {
            return new CoreMatchType(CoreMatchExtent.LEFT_PARTIAL, matchesRef);
        }

        if(matchesRight)
        {
            return new CoreMatchType(CoreMatchExtent.RIGHT_PARTIAL, matchesRef);
        }

        return null;
    }

    private boolean isLeftFlankMatch(final SAMRecord read)
    {
        String readString = read.getReadString();
        byte[] baseQuals = read.getBaseQualities();
        String leftFlankStr = mVariantReadContext.leftFlankStr();
        for(int i = 0; i < leftFlankStr.length(); ++i)
        {
            int readIndex = i + mVariantReadContext.AlignmentStart - read.getUnclippedStart();
            if(baseQuals[readIndex] < HIGH_BASE_QUAL_CUTOFF)
            {
                continue;
            }

            if(leftFlankStr.charAt(i) != readString.charAt(readIndex))
            {
                return false;
            }
        }

        return true;
    }

    private boolean isRightFlankMatch(final SAMRecord read)
    {
        String readString = read.getReadString();
        byte[] baseQuals = read.getBaseQualities();
        String rightFlankStr = mVariantReadContext.rightFlankStr();
        int rightFlankStartPos = mVariantReadContext.AlignmentStart + mVariantReadContext.CoreIndexEnd + 1;
        for(int i = 0; i < rightFlankStr.length(); ++i)
        {
            int readIndex = i + rightFlankStartPos - read.getUnclippedStart();
            if(baseQuals[readIndex] < HIGH_BASE_QUAL_CUTOFF)
            {
                continue;
            }

            if(rightFlankStr.charAt(i) != readString.charAt(readIndex))
            {
                return false;
            }
        }

        return true;
    }

}