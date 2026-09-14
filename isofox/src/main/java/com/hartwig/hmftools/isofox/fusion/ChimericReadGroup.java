package com.hartwig.hmftools.isofox.fusion;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.isofox.common.Read;

public class ChimericReadGroup
{
    private final List<Read> mReads;

    private boolean mComplete;

    public ChimericReadGroup(final Read read)
    {
        mReads = Lists.newArrayListWithCapacity(2);
        mReads.add(read);
        mComplete = readGroupComplete();
    }

    public ChimericReadGroup(final Read read1, final Read read2)
    {
        mReads = Lists.newArrayListWithCapacity(2);
        mReads.add(read1);
        mReads.add(read2);
        mComplete = readGroupComplete();
    }

    public final String id() { return mReads.get(0).Id; }

    public int size() { return mReads.size(); }

    public boolean isComplete() { return mComplete; }

    public List<Read> reads() { return mReads; }

    public void addRead(final Read read)
    {
        mReads.add(read);
        mComplete = readGroupComplete();
    }

    public boolean hasSuppAlignment() { return mReads.stream().anyMatch(x -> x.hasSuppAlignment()); }

    public String toString()
    {
        return String.format("%s reads(%d) complete(%s)", id(), mReads.size(), isComplete());
    }

    private boolean readGroupComplete()
    {
        int suppCount = 0;
        int nonSuppCount = 0;
        int expectedSuppCount = 0;
        int expectedNonSuppCount = 1;

        for(Read read : mReads)
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
