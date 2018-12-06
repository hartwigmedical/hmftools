package com.hartwig.hmftools.svanalysis.types;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.MIN_TEMPLATED_INSERTION_LENGTH;

import java.util.List;

public class SvLinkedPair {

    private SvVarData mFirst;
    private SvVarData mSecond;
    private boolean mFirstLinkOnStart;
    private boolean mSecondLinkOnStart;
    private String mLinkType;
    private int mLinkLength;
    private boolean mIsInferred;

    // other annotations
    private int mAssembledChainCount;
    private int mDBLenFirst;
    private int mDBLenSecond;
    private int mTraversedSVCount;
    private long mNearestSVDistance;
    private boolean mCopyNumberGain; // for TIs, is this from an additional fragment, not impacting the actual or derivative chromosomes
    private boolean mOnArmOfOrigin;

    public static final String LINK_TYPE_TI = "TI";
    public static final String LINK_TYPE_DB = "DB";
    public static final String LINK_TYPE_SGL = "SGL";

    public static String ASSEMBLY_MATCH_MATCHED = "MATCH";
    public static String ASSEMBLY_MATCH_DIFF = "DIFF";
    public static String ASSEMBLY_MATCH_INFER_ONLY = "INFER_ONLY";
    public static String ASSEMBLY_MATCH_NONE = "NONE";

    public SvLinkedPair(SvVarData first, SvVarData second, final String linkType, boolean firstLinkOnStart, boolean secondLinkOnStart)
    {
        mFirst = first;
        mSecond = second;
        mFirstLinkOnStart = firstLinkOnStart;
        mSecondLinkOnStart = secondLinkOnStart;
        mLinkType = linkType;
        mIsInferred = true;

        mAssembledChainCount = 0;
        mDBLenFirst = 0;
        mDBLenSecond = 0;
        mTraversedSVCount = 0;
        mNearestSVDistance = 0;
        mCopyNumberGain = false;
        mOnArmOfOrigin = false;

        int length = (int) (first.position(firstLinkOnStart) - second.position(secondLinkOnStart));
        mLinkLength = abs(length);

        if(mLinkType == LINK_TYPE_SGL)
        {
            if (mLinkLength < MIN_TEMPLATED_INSERTION_LENGTH)
            {
                mLinkLength = -mLinkLength;
            }
        }
        else
        {

            if (mLinkType == LINK_TYPE_TI && mLinkLength < MIN_TEMPLATED_INSERTION_LENGTH)
            {
                // re-label this as a deletion bridge
                mLinkType = LINK_TYPE_DB;
                mLinkLength = -mLinkLength;
            }
        }
    }

    public final SvVarData first() { return mFirst; }
    public final SvVarData second() { return mSecond; }
    public boolean firstLinkOnStart() { return mFirstLinkOnStart; }
    public boolean secondLinkOnStart() { return mSecondLinkOnStart; }
    public boolean firstUnlinkedOnStart() { return !mFirstLinkOnStart; }
    public boolean secondUnlinkedOnStart() { return !mSecondLinkOnStart; }

    public boolean getLinkedOnStart(final SvVarData var)
    {
        if(var.equals(mFirst))
            return mFirstLinkOnStart;
        else
            return mSecondLinkOnStart;
    }

    public final String linkType() { return mLinkType; }
    public final int length() { return mLinkLength; }

    public void setIsInferred(boolean toggle) { mIsInferred = toggle; }
    public boolean isInferred() { return mIsInferred; }
    public boolean isAssembled() { return !mIsInferred; }

    public void setAssembledChainCount(int count) { mAssembledChainCount = count; }
    public int getAssembledChainCount() { return mAssembledChainCount; }

    public void setDBLenFirst(int length) { mDBLenFirst = length; }
    public int getDBLenFirst() { return mDBLenFirst; }

    public void setDBLenSecond(int length) { mDBLenSecond = length; }
    public int getDBLenSecond() { return mDBLenSecond; }

    public void setNearestSVDistance(long distance) { mNearestSVDistance = distance; }
    public long getNearestSVDistance() { return mNearestSVDistance; }

    public void setCopyNumberGain(boolean isGain) { mCopyNumberGain = isGain; }
    public boolean hasCopyNumberGain() { return mCopyNumberGain; }

    public void setOnArmOfOrigin(boolean toggle) { mOnArmOfOrigin = toggle; }
    public boolean onArmOfOrigin() { return mOnArmOfOrigin; }

    public void setTraversedSVCount(int count) { mTraversedSVCount = count; }
    public int getTraversedSVCount() { return mTraversedSVCount; }

    public boolean hasVariantBE(final SvVarData var, boolean useStart)
    {
        return (var.equals(mFirst) && mFirstLinkOnStart == useStart || var.equals(mSecond) && mSecondLinkOnStart == useStart);
    }

    public boolean hasLinkClash(final SvLinkedPair otherPair)
    {
        return (hasVariantBE(otherPair.first(), otherPair.firstLinkOnStart())
            || hasVariantBE(otherPair.second(), otherPair.secondLinkOnStart()));
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
            if(pair.hasVariantBE(var, useStart))
                return pair;
        }

        return null;
    }

    public boolean matches(final SvLinkedPair other)
    {
        return this.matches(other, false);
    }

    public boolean matches(final SvLinkedPair other, boolean allowReplicated)
    {
        if(this == other)
            return true;

        // first and second can be in either order
        if(mFirst.equals(other.first(), allowReplicated) && mFirstLinkOnStart == other.firstLinkOnStart()
        && mSecond.equals(other.second(), allowReplicated) && mSecondLinkOnStart == other.secondLinkOnStart())
        {
            return true;
        }

        if(mFirst.equals(other.second(), allowReplicated) && mFirstLinkOnStart == other.secondLinkOnStart()
        && mSecond.equals(other.first(), allowReplicated) && mSecondLinkOnStart == other.firstLinkOnStart())
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
