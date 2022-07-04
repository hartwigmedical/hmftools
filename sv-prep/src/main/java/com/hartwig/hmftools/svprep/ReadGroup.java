package com.hartwig.hmftools.svprep;

import static java.lang.Math.min;

import java.util.List;

import com.google.common.collect.Lists;

public class ReadGroup
{
    private final List<ReadRecord> mReads;
    private int mMinPositionStart;
    private boolean mIsComplete;
    private boolean mHasSupplementary;

    public ReadGroup(final ReadRecord read)
    {
        mReads = Lists.newArrayListWithCapacity(2);
        addRead(read);
    }

    public List<ReadRecord> reads() { return mReads; }

    public void addRead(final ReadRecord read)
    {
        mReads.add(read);
        mHasSupplementary |= read.hasSuppAlignment() || read.isSupplementaryAlignment();

        if(mReads.size() == 1)
        {
            mMinPositionStart = read.start();
            mIsComplete = false;
            return;
        }

        mMinPositionStart = min(mMinPositionStart, read.start());
        mIsComplete = mReads.size() >= 3 || (mReads.size() == 2 && !mHasSupplementary);
    }

    public final String id() { return mReads.get(0).id(); }

    public int size() { return mReads.size(); }

    public int minStartPosition() { return mMinPositionStart; }

    public boolean isComplete() { return mIsComplete; }

    public boolean hasSuppAlignment() { return mHasSupplementary; }

    public String toString()
    {
        return String.format("%s reads(%d) minPos(%d) complete(%s)",
                id(), mReads.size(), mMinPositionStart, isComplete());
    }
}
