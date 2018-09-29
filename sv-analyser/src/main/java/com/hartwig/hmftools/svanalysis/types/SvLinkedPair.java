package com.hartwig.hmftools.svanalysis.types;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnalyser.MIN_TEMPLATED_INSERTION_LENGTH;

import java.util.List;

import com.google.common.collect.Lists;

public class SvLinkedPair {

    private SvClusterData mFirst;
    private SvClusterData mSecond;
    private boolean mFirstLinkOnStart;
    private boolean mSecondLinkOnStart;
    private String mLinkType;
    private int mLinkLength;
    private boolean mIsInferred;

    public static final String LINK_TYPE_TI = "TI";
    public static final String LINK_TYPE_DB = "DB";
    public static final String LINK_TYPE_SGL = "SGL";

    public SvLinkedPair(SvClusterData first, SvClusterData second, final String linkType, boolean firstLinkOnStart, boolean secondLinkOnStart)
    {
        mFirst = first;
        mSecond = second;
        mFirstLinkOnStart = firstLinkOnStart;
        mSecondLinkOnStart = secondLinkOnStart;
        mLinkType = linkType;
        mIsInferred = true;

        if(mLinkType == LINK_TYPE_SGL)
        {
            mLinkLength = (int)abs(first.position(true) - second.position(true));
        }
        else
        {
            int length = (int) (first.position(firstLinkOnStart) - second.position(secondLinkOnStart));
            mLinkLength = abs(length);

            if (mLinkType == LINK_TYPE_TI && mLinkLength < MIN_TEMPLATED_INSERTION_LENGTH)
            {
                // re-label this as a DB
                mLinkType = LINK_TYPE_DB;
                mLinkLength = -mLinkLength;
            }
        }
    }

    public final SvClusterData first() { return mFirst; }
    public final SvClusterData second() { return mSecond; }
    public boolean firstLinkOnStart() { return mFirstLinkOnStart; }
    public boolean secondLinkOnStart() { return mSecondLinkOnStart; }
    public boolean firstUnlinkedOnStart() { return !mFirstLinkOnStart; }
    public boolean secondUnlinkedOnStart() { return !mSecondLinkOnStart; }

    public final String linkType() { return mLinkType; }
    public final int length() { return mLinkLength; }

    public void setIsInferred(boolean toggle) { mIsInferred = toggle; }
    public boolean isInferred() { return mIsInferred; }

    public boolean hasVariantBE(final SvClusterData var, boolean useStart)
    {
        return (var.equals(mFirst) && mFirstLinkOnStart == useStart || var.equals(mSecond) && mSecondLinkOnStart == useStart);
    }

    public boolean hasLinkClash(final SvLinkedPair otherPair)
    {
        return (hasVariantBE(otherPair.first(), otherPair.firstLinkOnStart())
            || hasVariantBE(otherPair.second(), otherPair.secondLinkOnStart()));
    }

    public final String toString()
    {
        boolean firstLinkBE = mLinkType == LINK_TYPE_SGL ? !mFirstLinkOnStart : mFirstLinkOnStart;
        boolean secondLinkBE = mLinkType == LINK_TYPE_SGL ? !mSecondLinkOnStart : mSecondLinkOnStart;

        return String.format("%s %s:%d:%s & %s %s:%d:%s",
                first().id(), first().chromosome(mFirstLinkOnStart), first().position(firstLinkBE), firstLinkBE ? "start":"end",
                second().id(), second().chromosome(mSecondLinkOnStart), second().position(secondLinkBE), secondLinkBE ? "start":"end");

    }

    public static final SvLinkedPair findLinkedPair(final List<SvLinkedPair> linkedPairs, final SvClusterData var, boolean useStart)
    {
        for(final SvLinkedPair pair : linkedPairs)
        {
            if(var.equals(pair.first()) && useStart == pair.firstLinkOnStart())
                return pair;

            if(var.equals(pair.second()) && useStart == pair.secondLinkOnStart())
                return pair;
        }

        return null;
    }

    public static String ASSEMBLY_MATCH_MATCHED = "MATCH";
    public static String ASSEMBLY_MATCH_DIFF = "DIFF";
    public static String ASSEMBLY_MATCH_ASMB_ONLY = "ASMB_ONLY";
    public static String ASSEMBLY_MATCH_LINK_ONLY = "LINK_ONLY";
    public static String ASSEMBLY_MATCH_NONE = "NONE";

    public String getAssemblyMatchType(final SvClusterData var)
    {
        final String firstAssembly = mFirstLinkOnStart ? mFirst.getAssemblyStart() : mFirst.getAssemblyEnd();
        final String secondAssembly = mSecondLinkOnStart ? mSecond.getAssemblyStart() : mSecond.getAssemblyEnd();

        if((var == mFirst && firstAssembly.isEmpty()) || (var == mSecond && secondAssembly.isEmpty()))
        {
            return ASSEMBLY_MATCH_LINK_ONLY;
        }

        if(firstAssembly.equals(secondAssembly))
        {
            return ASSEMBLY_MATCH_MATCHED;
        }

        String[] firstAssemblyList = firstAssembly.split(";");
        String[] secondAssemblyList = secondAssembly.split(";");

        for(int i = 0; i < firstAssemblyList.length; ++i)
        {
            for(int j = 0; j < secondAssemblyList.length; ++j)
            {
                if(firstAssemblyList[i].equals(secondAssemblyList[j]))
                    return ASSEMBLY_MATCH_MATCHED;
            }
        }


        return ASSEMBLY_MATCH_DIFF;
    }

}
