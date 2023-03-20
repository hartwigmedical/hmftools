package com.hartwig.hmftools.isofox.fusion;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.isofox.common.ReadRecord;

public class ChimericReadGroup
{
    private final List<ReadRecord> mReads;

    private boolean mComplete;

    public ChimericReadGroup(final ReadRecord read)
    {
        mReads = Lists.newArrayListWithCapacity(2);
        mReads.add(read);
        mComplete = readGroupComplete(mReads);
    }

    public ChimericReadGroup(final ReadRecord read1, final ReadRecord read2)
    {
        mReads = Lists.newArrayListWithCapacity(2);
        mReads.add(read1);
        mReads.add(read2);
        mComplete = readGroupComplete(mReads);
    }

    public final String id() { return mReads.get(0).Id; }

    public int size() { return mReads.size(); }

    public boolean isComplete() { return mComplete; }

    public List<ReadRecord> reads() { return mReads; }

    public void addRead(final ReadRecord read)
    {
        mReads.add(read);
        mComplete = readGroupComplete(mReads);
    }

    public boolean hasSuppAlignment() { return mReads.stream().anyMatch(x -> x.hasSuppAlignment()); }

    public String toString()
    {
        return String.format("%s reads(%d) complete(%s)", id(), mReads.size(), isComplete());
    }

    public static boolean readGroupComplete(final List<ReadRecord> reads)
    {
        int suppCount = 0;
        int nonSuppCount = 0;
        int expectedSuppCount = 0;
        int expectedNonSuppCount = 1;

        for(ReadRecord read : reads)
        {
            if(read.isReadPaired() && !read.isMateUnmapped())
            {
                expectedNonSuppCount = 2;
            }

            if(read.isSupplementaryAlignment())
            {
                ++suppCount;
            }
            else
            {
                ++nonSuppCount;

                if(read.hasSuppAlignment())
                {
                    ++expectedSuppCount;
                }
            }
        }

        return (expectedNonSuppCount == nonSuppCount) && (expectedSuppCount == suppCount);
    }
}
