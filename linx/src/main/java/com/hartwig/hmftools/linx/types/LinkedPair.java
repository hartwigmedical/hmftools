package com.hartwig.hmftools.linx.types;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;

import java.util.List;

public class LinkedPair
{

    private SvBreakend mFirstBreakend;
    private SvBreakend mSecondBreakend;
    private int mLinkLength;
    private boolean mIsInferred;

    // other annotations
    private String mLinkReason;
    private int mLinkIndex;
    private int mTraversedSVCount;
    private int mNextSvDistance;
    private int mNextClusteredSvDistance;
    private String mLocationType;
    private int mOverlapCount;
    private boolean mHasCopyNumberGain;
    private int mIndelCount;

    private String mExonMatchData;

    public static final String LOCATION_TYPE_UNCLEAR = "Unclear";
    public static final String LOCATION_TYPE_REMOTE = "Remote"; // TI is not on arm with any chain end
    public static final String LOCATION_TYPE_EXTERNAL = "External"; // TI is on arm with a chain end but outside its bounds
    public static final String LOCATION_TYPE_INTERNAL = "Internal";

    public LinkedPair(final SvBreakend first, final SvBreakend second)
    {
        mFirstBreakend = first;
        mSecondBreakend = second;
        mIsInferred = true;

        mLinkReason = "";
        mLinkIndex = -1;
        mTraversedSVCount = 0;
        mNextSvDistance = 0;
        mNextClusteredSvDistance = 0;
        mLocationType = LOCATION_TYPE_UNCLEAR;
        mOverlapCount = 0;
        mHasCopyNumberGain = false;
        mIndelCount = 0;
        mExonMatchData = "";

        int length = first.position() - second.position();
        mLinkLength = abs(length);
    }

    public static LinkedPair from(SvVarData first, SvVarData second, boolean firstLinkOnStart, boolean secondLinkOnStart)
    {
        return new LinkedPair(first.getBreakend(firstLinkOnStart), second.getBreakend(secondLinkOnStart));
    }

    public static LinkedPair from(final SvBreakend first, final SvBreakend second)
    {
        return new LinkedPair(first, second);
    }

    public static LinkedPair copy(final LinkedPair other)
    {
        LinkedPair newPair = LinkedPair.from(other.firstBreakend(), other.secondBreakend());
        newPair.setLinkReason(other.getLinkReason(), other.getLinkIndex());

        if(other.isAssembled())
            newPair.setIsAssembled();

        return newPair;
    }

    public final SvVarData first() { return mFirstBreakend.getSV(); }
    public final SvVarData second() { return mSecondBreakend.getSV(); }

    public boolean firstLinkOnStart() { return mFirstBreakend.usesStart(); }
    public boolean secondLinkOnStart() { return mSecondBreakend.usesStart(); }
    public boolean firstUnlinkedOnStart() { return !mFirstBreakend.usesStart(); }
    public boolean secondUnlinkedOnStart() { return !mSecondBreakend.usesStart(); }

    public final SvBreakend getBreakend(int se) { return getBreakend(isStart(se)); }

    public final SvBreakend getBreakend(boolean isStart)
    {
        // finds the earlier breakend of the 2, ie with the lower position
        // unless explicitly linked on a SGL mapping, INFs and SGLs return their only breakend (ie the first)
        final SvBreakend beFirst = mFirstBreakend.getSV().isSglBreakend() && mFirstBreakend.usesStart() ?
                mFirstBreakend.getSV().getBreakend(true) : mFirstBreakend;

        final SvBreakend beSecond = mSecondBreakend.getSV().isSglBreakend() && mSecondBreakend.usesStart() ?
                mSecondBreakend.getSV().getBreakend(true) : mSecondBreakend;

        if(isStart)
            return beFirst.position() <= beSecond.position() ? beFirst : beSecond;
        else
            return beFirst.position() > beSecond.position() ? beFirst : beSecond;
    }

    public final SvBreakend firstBreakend() { return mFirstBreakend; }
    public final SvBreakend secondBreakend() { return mSecondBreakend; }

    public final String chromosome() { return mFirstBreakend.chromosome(); }

    public final int positionDistance() { return mLinkLength; }
    public final int baseLength() { return mLinkLength + 1; }

    public void setIsAssembled() { mIsInferred = false; }
    public boolean isInferred() { return mIsInferred; }
    public boolean isAssembled() { return !mIsInferred; }

    public String getLinkReason() { return mLinkReason; }
    public void setLinkReason(String reason, int index)
    {
        mLinkReason = reason;
        mLinkIndex = index;
    }

    public int getLinkIndex() { return mLinkIndex; }

    public void setNextSVData(int distance, int clusterDistance)
    {
        mNextSvDistance = distance;
        mNextClusteredSvDistance = clusterDistance;
    }

    public int getNextSvDistance() { return mNextSvDistance; }
    public int getNextClusteredSvDistance() { return mNextClusteredSvDistance; }

    public void setLocationType(final String type) { mLocationType = type; }
    public final String locationType() { return mLocationType; }

    public void setOverlapCount(int count) { mOverlapCount = count; }
    public int overlapCount() { return mOverlapCount; }

    public void setHasCopyNumberGain(boolean toggle) { mHasCopyNumberGain = toggle; }
    public boolean hasCopyNumberGain() { return mHasCopyNumberGain; }

    public void setTraversedSVCount(int count) { mTraversedSVCount = count; }
    public int getTraversedSVCount() { return mTraversedSVCount; }

    public void setExonMatchData(final String data) { mExonMatchData = data; }
    public final String getExonMatchData() { return mExonMatchData; }

    public boolean hasBreakend(final SvVarData var, boolean useStart)
    {
        return (var == mFirstBreakend.getSV() && mFirstBreakend.usesStart() == useStart)
                || (var == mSecondBreakend.getSV() && mSecondBreakend.usesStart() == useStart);
    }

    public boolean hasBreakend(final SvBreakend breakend)
    {
        return hasBreakend(breakend.getSV(), breakend.usesStart());
    }

    public boolean hasLinkClash(final LinkedPair otherPair)
    {
        return (hasBreakend(otherPair.first(), otherPair.firstLinkOnStart())
            || hasBreakend(otherPair.second(), otherPair.secondLinkOnStart()));
    }

    public void switchSVs()
    {
        final SvBreakend tmp = mSecondBreakend;
        mSecondBreakend = mFirstBreakend;
        mFirstBreakend = tmp;
    }

    public final String toString()
    {
        return String.format("%s %s:%d:%s & %s %s:%d:%s",
                mFirstBreakend.getSV().id(), mFirstBreakend.chromosome(), mFirstBreakend.position(),
                mFirstBreakend.usesStart() ? "start" : "end",
                mSecondBreakend.getSV().id(), mSecondBreakend.chromosome(), mSecondBreakend.position(),
                mSecondBreakend.usesStart() ? "start" : "end");
    }

    public boolean matches(final LinkedPair other)
    {
        if(this == other)
            return true;

        // first and second can be in either order
        if(mFirstBreakend == other.firstBreakend() && mSecondBreakend == other.secondBreakend())
            return true;

        if(mFirstBreakend == other.secondBreakend() && mSecondBreakend == other.firstBreakend())
            return true;

        return false;
    }

    public boolean oppositeMatch(final LinkedPair other)
    {
        if(this == other)
            return true;

        if(mFirstBreakend.getSV() == other.first() && mFirstBreakend.usesStart() != other.firstLinkOnStart()
        && mSecondBreakend.getSV() ==other.second() && mSecondBreakend.usesStart() != other.secondLinkOnStart())
        {
            return true;
        }

        if(mFirstBreakend.getSV() == other.second() && mFirstBreakend.usesStart() != other.secondLinkOnStart()
        && mSecondBreakend.getSV() == other.first() && mSecondBreakend.usesStart() != other.firstLinkOnStart())
        {
            return true;
        }

        return false;
    }

    public boolean sameVariants(final LinkedPair other)
    {
        // first and second can be in either order
        if(mFirstBreakend.getSV() == other.first() && mSecondBreakend.getSV() == other.second())
            return true;

        if(mFirstBreakend.getSV() == other.second() && mSecondBreakend.getSV() == other.first())
            return true;

        return false;
    }

    public final SvVarData getOtherSV(final SvVarData var) { return mFirstBreakend.getSV() == var ? mSecondBreakend.getSV() : mFirstBreakend.getSV(); }

    public final SvBreakend getOtherBreakend(final SvBreakend breakend)
    {
        final SvBreakend lower = getBreakend(true);
        return breakend == lower ? getBreakend(false) : lower;
    }

    public boolean hasVariant(final SvVarData var) { return mFirstBreakend.getSV() == var || mSecondBreakend.getSV() == var; }

    public boolean isDupLink() { return mFirstBreakend.getSV() == mSecondBreakend.getSV() && mFirstBreakend.type() == DUP; }

    public static boolean hasLinkClash(final List<LinkedPair> links1, final List<LinkedPair> links2)
    {
        for(final LinkedPair pair : links1)
        {
            if(links2.stream().anyMatch(x -> pair.hasLinkClash(x)))
                return true;
        }

        return false;
    }

}
