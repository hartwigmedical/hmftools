package com.hartwig.hmftools.bamtools.metrics;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;

import java.util.List;

import com.google.common.collect.Lists;

import htsjdk.samtools.SAMRecord;

public class ReadGroup
{
    public final boolean IsConsensus;
    private final List<SAMRecord> mReads;
    private boolean mAllReadsPresent;

    public ReadGroup(final SAMRecord read, boolean isConsensus)
    {
        IsConsensus = isConsensus;
        mReads = Lists.newArrayListWithCapacity(2);
        addRead(read);
        mAllReadsPresent = false;
        checkComplete();
    }

    public final String id() { return mReads.get(0).getReadName(); }
    public List<SAMRecord> reads() { return mReads; }
    public int size() { return mReads.size(); }
    public boolean allReadsPresent() { return mAllReadsPresent; }

    public void addRead(final SAMRecord read)
    {
        mReads.add(read);
        checkComplete();
    }

    private void checkComplete()
    {
        if(mAllReadsPresent)
            return;

        int suppCount = 0;
        int nonSuppCount = 0;
        int expectedSuppCount = 0;
        int expectedNonSuppCount = 1;

        for(SAMRecord read : mReads)
        {
            if(read.getReadPairedFlag() && !read.getMateUnmappedFlag())
            {
                expectedNonSuppCount = 2;
            }

            if(read.getSupplementaryAlignmentFlag())
            {
                ++suppCount;
            }
            else
            {
                ++nonSuppCount;
            }

            if(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
            {
                if(!read.getSupplementaryAlignmentFlag())
                {
                    ++expectedSuppCount;
                }
            }
        }

        mAllReadsPresent = (expectedNonSuppCount == nonSuppCount) && (expectedSuppCount == suppCount);
    }

    public String toString()
    {
        return String.format("reads(%d) initReadStart(%s:%d) id(%s)",
                mReads.size(), mReads.get(0).getContig(), mReads.get(0).getAlignmentStart(), id());
    }
}
