package com.hartwig.hmftools.svanalysis.types;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.MIN_TEMPLATED_INSERTION_LENGTH;

import java.util.List;

public class SvLinkedPair {

    private int mId;
    private SvVarData mFirst;
    private SvVarData mSecond;
    private boolean mFirstLinkOnStart;
    private boolean mSecondLinkOnStart;
    private String mLinkType;
    private int mLinkLength;
    private boolean mIsInferred;

    // other annotations
    private int mDBLenFirst;
    private int mDBLenSecond;
    private int mTraversedSVCount;
    private int mNextSVDistance;
    private int mNextSVTraversedCount;
    private boolean mCopyNumberGain; // for TIs, is this from an additional fragment, not impacting the actual or derivative chromosomes
    private boolean mOnArmOfOrigin;

    private String mExonMatchData;

    public static final String LINK_TYPE_TI = "TI";
    public static final String LINK_TYPE_DB = "DB";

    public static String ASSEMBLY_MATCH_MATCHED = "MATCH";
    public static String ASSEMBLY_MATCH_DIFF = "DIFF";
    public static String ASSEMBLY_MATCH_INFER_ONLY = "INFER_ONLY";
    public static String ASSEMBLY_MATCH_NONE = "NONE";

    public SvLinkedPair(SvVarData first, SvVarData second, final String linkType, boolean firstLinkOnStart, boolean secondLinkOnStart)
    {
        mId = 0;
        mFirst = first;
        mSecond = second;
        mFirstLinkOnStart = firstLinkOnStart;
        mSecondLinkOnStart = secondLinkOnStart;
        mLinkType = linkType;
        mIsInferred = true;

        mDBLenFirst = 0;
        mDBLenSecond = 0;
        mTraversedSVCount = 0;
        mNextSVDistance = 0;
        mNextSVTraversedCount = 0;
        mCopyNumberGain = false;
        mOnArmOfOrigin = false;
        mExonMatchData = "";

        int length = (int) (first.position(firstLinkOnStart) - second.position(secondLinkOnStart));
        mLinkLength = abs(length);

        if (mLinkType == LINK_TYPE_TI && mLinkLength < MIN_TEMPLATED_INSERTION_LENGTH)
        {
            // re-label this as a deletion bridge
            mLinkType = LINK_TYPE_DB;
            mLinkLength = -mLinkLength;
        }

        // adjust the length of DBs to reflect the position convention for opposite breakend orientations
        if(mLinkType == LINK_TYPE_DB)
            --mLinkLength;
    }

    public static SvLinkedPair from(final SvBreakend first, final SvBreakend second, final String linkType)
    {
        return new SvLinkedPair(first.getSV(), second.getSV(), linkType, first.usesStart(), second.usesStart());
    }


    public int id() { return mId; }
    public void setId(int Id) { mId = 0; }

    public final SvVarData first() { return mFirst; }
    public final SvVarData second() { return mSecond; }

    public void replaceFirst(final SvVarData var) { mFirst = var; }
    public void replaceSecond(final SvVarData var) { mSecond = var; }

    public boolean firstLinkOnStart() { return mFirstLinkOnStart; }
    public boolean secondLinkOnStart() { return mSecondLinkOnStart; }
    public boolean firstUnlinkedOnStart() { return !mFirstLinkOnStart; }
    public boolean secondUnlinkedOnStart() { return !mSecondLinkOnStart; }

    public final SvBreakend getBreakend(boolean isStart)
    {
        // finds the earlier breakend of the 2, ie with the lower position
        final SvBreakend beFirst = mFirst.type() == SGL ? mFirst.getBreakend(true) : mFirst.getBreakend(firstLinkOnStart());
        final SvBreakend beSecond = mSecond.type() == SGL ? mSecond.getBreakend(true) : mSecond.getBreakend(secondLinkOnStart());

        if(isStart)
            return beFirst.position() < beSecond.position() ? beFirst : beSecond;
        else
            return beFirst.position() > beSecond.position() ? beFirst : beSecond;
    }

    public final String linkType() { return mLinkType; }
    public final int length() { return mLinkLength; }

    public void setIsInferred(boolean toggle) { mIsInferred = toggle; }
    public boolean isInferred() { return mIsInferred; }
    public boolean isAssembled() { return !mIsInferred; }
    public String assemblyInferredStr() { return mIsInferred ? "inferred" : "assembly"; }

    public void setDBLenFirst(int length) { mDBLenFirst = length; }
    public int getDBLenFirst() { return mDBLenFirst; }

    public void setDBLenSecond(int length) { mDBLenSecond = length; }
    public int getDBLenSecond() { return mDBLenSecond; }

    public void setNextSVData(int distance, int traversedCount)
    {
        mNextSVDistance = distance;
        mNextSVTraversedCount = traversedCount;
    }

    public long getNextSVDistance() { return mNextSVDistance; }
    public long getNextSVTraversedCount() { return mNextSVTraversedCount; }

    public void setCopyNumberGain(boolean isGain) { mCopyNumberGain = isGain; }
    public boolean hasCopyNumberGain() { return mCopyNumberGain; }

    public void setOnArmOfOrigin(boolean toggle) { mOnArmOfOrigin = toggle; }
    public boolean onArmOfOrigin() { return mOnArmOfOrigin; }

    public void setTraversedSVCount(int count) { mTraversedSVCount = count; }
    public int getTraversedSVCount() { return mTraversedSVCount; }

    public void setExonMatchData(final String data) { mExonMatchData = data; }
    public final String getExonMatchData() { return mExonMatchData; }


    public boolean hasBreakend(final SvVarData var, boolean useStart, boolean allowReplicated)
    {
        return (var.equals(mFirst, allowReplicated) && mFirstLinkOnStart == useStart
                || var.equals(mSecond, allowReplicated) && mSecondLinkOnStart == useStart);
    }

    public boolean hasBreakend(final SvVarData var, boolean useStart)
    {
        return hasBreakend(var, useStart, false);
    }

    public boolean hasBreakend(final SvBreakend breakend, boolean allowReplicated)
    {
        return hasBreakend(breakend.getSV(), breakend.usesStart(), allowReplicated);
    }

    public boolean hasBreakend(final SvBreakend breakend)
    {
        return hasBreakend(breakend.getSV(), breakend.usesStart());
    }

    public boolean hasLinkClash(final SvLinkedPair otherPair)
    {
        return (hasBreakend(otherPair.first(), otherPair.firstLinkOnStart())
            || hasBreakend(otherPair.second(), otherPair.secondLinkOnStart()));
    }

    public void switchSVs()
    {
        final SvVarData tmp = mSecond;
        mSecond = mFirst;
        mFirst = tmp;

        boolean tmp2 = mSecondLinkOnStart;
        mSecondLinkOnStart = mFirstLinkOnStart;
        mFirstLinkOnStart = tmp2;
    }

    public final String toString()
    {
        return svToString(mFirst, mFirstLinkOnStart) + " & " + svToString(mSecond, mSecondLinkOnStart);
    }

    private static final String svToString(final SvVarData var, boolean linkedOnStart)
    {
        if(var.type() != SGL)
        {
            return String.format("%s %s:%d:%s",
                    var.id(), var.chromosome(linkedOnStart), var.position(linkedOnStart), linkedOnStart ? "start" : "end");
        }

        if(linkedOnStart)
        {
            return String.format("%s %s:%d SGL-on-known",
                    var.id(), var.chromosome(true), var.position(true));
        }
        else
        {
            return String.format("%s %s:%d SGL-on-null",
                    var.id(), var.chromosome(true), var.position(true));
        }


    }

    public static final SvLinkedPair findLinkedPair(final List<SvLinkedPair> linkedPairs, final SvVarData var, boolean useStart)
    {
        for(final SvLinkedPair pair : linkedPairs)
        {
            if(pair.hasBreakend(var, useStart))
                return pair;
        }

        return null;
    }

    public boolean matches(final SvLinkedPair other)
    {
        if(this == other)
            return true;

        // first and second can be in either order
        if(mFirst.equals(other.first(), true) && mFirstLinkOnStart == other.firstLinkOnStart()
        && mSecond.equals(other.second(), true) && mSecondLinkOnStart == other.secondLinkOnStart())
        {
            return true;
        }

        if(mFirst.equals(other.second(), true) && mFirstLinkOnStart == other.secondLinkOnStart()
        && mSecond.equals(other.first(), true) && mSecondLinkOnStart == other.firstLinkOnStart())
        {
            return true;
        }

        return false;
    }

    public boolean oppositeMatch(final SvLinkedPair other)
    {
        if(this == other)
            return true;

        if(mFirst.equals(other.first(), true) && mFirstLinkOnStart != other.firstLinkOnStart()
        && mSecond.equals(other.second(), true) && mSecondLinkOnStart != other.secondLinkOnStart())
        {
            return true;
        }

        if(mFirst.equals(other.second(), true) && mFirstLinkOnStart != other.secondLinkOnStart()
        && mSecond.equals(other.first(), true) && mSecondLinkOnStart != other.firstLinkOnStart())
        {
            return true;
        }

        return false;
    }

    public boolean sameVariants(final SvLinkedPair other)
    {
        // first and second can be in either order
        if(mFirst.equals(other.first(), true) && mSecond.equals(other.second(), true))
            return true;

        if(mFirst.equals(other.second(), true) && mSecond.equals(other.first(), true))
            return true;

        return false;
    }

    public final SvVarData getOtherSV(final SvVarData var)
    {
        return mFirst.equals(var, true) ? mSecond : mFirst;
    }

    public final SvBreakend getOtherBreakend(final SvBreakend breakend)
    {
        final SvBreakend lower = getBreakend(true);
        return breakend == lower ? getBreakend(false) : lower;
    }

    public boolean hasAnySameVariant(final SvLinkedPair other)
    {
        return (mFirst.equals(other.first(), true) || mSecond.equals(other.second(), true)
            || mFirst.equals(other.second(), true) || mSecond.equals(other.first(), true));
    }

    public static void removedLinksWithSV(List<SvLinkedPair> list, final SvVarData var)
    {
        int index = 0;
        while(index < list.size())
        {
            SvLinkedPair pair = list.get(index);

            if(pair.first() == var || pair.second() == var)
            {
                list.remove(index);
            }
            else
            {
                ++index;
            }
        }
    }


}
