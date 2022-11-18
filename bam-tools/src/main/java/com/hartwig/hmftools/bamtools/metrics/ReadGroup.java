package com.hartwig.hmftools.bamtools.metrics;

import java.util.List;

import com.google.common.collect.Lists;

import htsjdk.samtools.SAMRecord;

public class ReadGroup
{
    private final List<SAMRecord> mReads;

    public ReadGroup(final SAMRecord read)
    {
        mReads = Lists.newArrayListWithCapacity(2);
        addRead(read);
    }

    public final String id() { return mReads.get(0).getReadName(); }
    public List<SAMRecord> reads() { return mReads; }
    public int size() { return mReads.size(); }

    public void addRead(final SAMRecord read)
    {
        mReads.add(read);
    }

    public String toString()
    {
        return String.format("reads(%d) initReadStart(%s:%d) id(%s)",
                mReads.size(), mReads.get(0).getContig(), mReads.get(0).getAlignmentStart(), id());
    }
}
