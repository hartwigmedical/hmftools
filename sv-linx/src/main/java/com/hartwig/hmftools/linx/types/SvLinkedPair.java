package com.hartwig.hmftools.linx.types;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.chaining.LinkFinder.getMinTemplatedInsertionLength;
import static com.hartwig.hmftools.linx.types.LinkType.DELETION_BRIDGE;
import static com.hartwig.hmftools.linx.types.LinkType.TEMPLATED_INSERTION;
import static com.hartwig.hmftools.linx.types.LinkType.linkTypeStr;

import java.util.List;

import com.hartwig.hmftools.linx.visualiser.data.Link;

public class SvLinkedPair {

    private SvBreakend mFirstBreakend;
    private SvBreakend mSecondBreakend;
    private LinkType mLinkType;
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

    public static String LOCATION_TYPE_UNCLEAR = "Unclear";
    public static String LOCATION_TYPE_REMOTE = "Remote"; // TI is not on arm with any chain end
    public static String LOCATION_TYPE_EXTERNAL = "External"; // TI is on arm with a chain end but outside its bounds
    public static String LOCATION_TYPE_INTERNAL = "Internal";

    public SvLinkedPair(final SvBreakend first, final SvBreakend second, final LinkType linkType)
    {
        mFirstBreakend = first;
        mSecondBreakend = second;
        mLinkType = linkType;
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

        int minTILength = getMinTemplatedInsertionLength(first, second);

        if (mLinkType == TEMPLATED_INSERTION && mLinkLength < minTILength)
        {
            // re-label this as a deletion bridge and give it a negative length to show the overlap
            mLinkType = DELETION_BRIDGE;
            mLinkLength = -mLinkLength;

            // if the link is marked as assembled later on, this is reversed
        }

        // adjust the length of DBs to reflect the position convention for opposite breakend orientations
        if(mLinkType == DELETION_BRIDGE)
            --mLinkLength;
    }

    public static SvLinkedPair from(SvVarData first, SvVarData second, final LinkType linkType, boolean firstLinkOnStart, boolean secondLinkOnStart)
    {
        return new SvLinkedPair(first.getBreakend(firstLinkOnStart), second.getBreakend(secondLinkOnStart), linkType);
    }

    public static SvLinkedPair from(final SvBreakend first, final SvBreakend second)
    {
        return new SvLinkedPair(first, second, TEMPLATED_INSERTION);
    }

    public static SvLinkedPair copy(final SvLinkedPair other)
    {
        SvLinkedPair newPair = SvLinkedPair.from(other.firstBreakend(), other.secondBreakend());
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

    public final LinkType linkType() { return mLinkType; }
    public final int length() { return mLinkLength; }

    public void setIsAssembled()
    {
        mIsInferred = false;

        if(mLinkLength < 0)
        {
            mLinkType = TEMPLATED_INSERTION;
            mLinkLength = -mLinkLength;
        }
    }
    public boolean isInferred() { return mIsInferred; }
    public boolean isAssembled() { return !mIsInferred; }
    public String assemblyInferredStr() { return mIsInferred ? "Inferred" : "Assembly"; }

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

    public void setIndelCount(final int count) { mIndelCount = count; }
    public int getIndelCount() { return mIndelCount; }

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

    public boolean hasLinkClash(final SvLinkedPair otherPair)
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

    private static final String svToString(final SvVarData var, boolean linkedOnStart)
    {
        if(!var.isSglBreakend())
        {
            return String.format("%s %s:%d:%s",
                    var.id(), var.chrShort(linkedOnStart), var.position(linkedOnStart),
                    linkedOnStart ? "start" : "end");
        }
        else
        {
            return String.format("%s %s:%d SGL",
                    var.id(), var.chrShort(true), var.position(true));
        }
    }

    public boolean matches(final SvLinkedPair other)
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

    public boolean oppositeMatch(final SvLinkedPair other)
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

    public boolean sameVariants(final SvLinkedPair other)
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

    public boolean isDupLink() { return mFirstBreakend.getSV() == mSecondBreakend.getSV() && mFirstBreakend.getSV().type() == DUP; }

    public static boolean hasLinkClash(final List<SvLinkedPair> links1, final List<SvLinkedPair> links2)
    {
        for(final SvLinkedPair pair : links1)
        {
            if(links2.stream().anyMatch(x -> pair.hasLinkClash(x)))
                return true;
        }

        return false;
    }

}
