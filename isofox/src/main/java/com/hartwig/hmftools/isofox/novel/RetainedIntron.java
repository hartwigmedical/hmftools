package com.hartwig.hmftools.isofox.novel;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.Strings.appendStrList;

import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.isofox.common.RegionReadData;
import com.hartwig.hmftools.isofox.common.TransExonRef;

public class RetainedIntron
{
    private final List<RegionReadData> mRegions;
    private final boolean mIsStart;
    private final int mPosition;

    private int mFragmentCount;
    private int mSplicedFragmentCount;
    private int mReadDepth;

    private final String mTranscriptInfo;

    public RetainedIntron(final List<RegionReadData> regions, boolean isStart)
    {
        mRegions = regions;
        mIsStart = isStart;

        mFragmentCount = 0;
        mSplicedFragmentCount = 0;
        mReadDepth = 0;

        mPosition = isStart ? regions.get(0).start() : regions.get(0).end();

        StringJoiner trancriptInfo = new StringJoiner(";");

        for(RegionReadData region : mRegions)
        {
            for (final TransExonRef transRef : region.getTransExonRefs())
            {
                trancriptInfo.add(format("%s-%s", transRef.TransName, transRef.ExonRank));
            }
        }

        mTranscriptInfo = trancriptInfo.toString();
    }

    public final List<RegionReadData> regions() { return mRegions; }
    public boolean isStart() { return mIsStart; }
    public int position() { return mPosition; }

    public String type(boolean isForwardStrand) { return (isForwardStrand == mIsStart) ? "ACCEPTOR" : "DONOR"; }

    public String transcriptInfo() { return mTranscriptInfo; }

    public void addFragmentCount(boolean spliceSupport)
    {
        ++mFragmentCount;

        if(spliceSupport)
            ++mSplicedFragmentCount;
    }

    public int getFragmentCount() { return mFragmentCount; }
    public int getSplicedFragmentCount() { return mSplicedFragmentCount; }

    public void addReadDepth() { ++mReadDepth; }
    public void setReadDepth(int count) { mReadDepth = count; }
    public int getDepth() { return mReadDepth; }

    public boolean matches(final RetainedIntron other)
    {
        return other.position() == mPosition && other.isStart() == mIsStart;
    }

    public String toString()
    {
        return format("pos(%s @ %s) regions(%d)", mPosition, mIsStart ? "start" : "end", mRegions.size());
    }
}
