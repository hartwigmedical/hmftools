package com.hartwig.hmftools.svanalysis.types;

import static java.lang.Math.abs;

public class SvLinkedPair {

    private SvClusterData mFirst;
    private SvClusterData mSecond;
    private boolean mFirstLinkOnStart;
    private boolean mSecondLinkOnStart;
    private String mLinkType;
    private int mLinkLength;

    public static final String LINK_TYPE_TI = "TI";
    public static final String LINK_TYPE_DB = "DB";
    public static final String LINK_TYPE_DUP_BE = "DBE";

    public SvLinkedPair(SvClusterData first, SvClusterData second, final String linkType, boolean firstLinkOnStart, boolean secondLinkOnStart)
    {
        mFirst = first;
        mSecond = second;
        mFirstLinkOnStart = firstLinkOnStart;
        mSecondLinkOnStart = secondLinkOnStart;
        mLinkType = linkType;

        int length = (int)(first.position(firstLinkOnStart) - second.position(secondLinkOnStart));
        mLinkLength = abs(length);
    }

    public final SvClusterData first() { return mFirst; }
    public final SvClusterData second() { return mSecond; }
    public boolean firstLinkOnStart() { return mFirstLinkOnStart; }
    public boolean secondLinkOnStart() { return mSecondLinkOnStart; }
    public final String linkType() { return mLinkType; }
    public final int length() { return mLinkLength; }

    public boolean hasVariantBE(final SvClusterData var, boolean useStart)
    {
        return (var.equals(mFirst) && mFirstLinkOnStart == useStart || var.equals(mSecond) && mSecondLinkOnStart == useStart);
    }

    public final String toString()
    {
        return String.format("%s %s:%d:%s & %s %s:%d:%s",
                first().id(), first().chromosome(mFirstLinkOnStart), first().position(mFirstLinkOnStart), mFirstLinkOnStart ? "start":"end",
                second().id(), second().chromosome(mSecondLinkOnStart), second().position(mSecondLinkOnStart), mSecondLinkOnStart ? "start":"end");
    }

}
