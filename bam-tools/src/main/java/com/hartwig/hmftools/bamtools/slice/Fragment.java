package com.hartwig.hmftools.bamtools.slice;

import static java.lang.String.format;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.region.BasePosition;

import htsjdk.samtools.SAMRecord;

public class Fragment
{
    private final String mReadId;
    private final boolean mUnpaired;

    private final List<ReadInfo> mReceivedReads;
    private final List<ReadInfo> mPendingReads;

    private final List<SAMRecord> mCachedReads;

    public Fragment(final SAMRecord read)
    {
        mReadId = read.getReadName();
        mReceivedReads = Lists.newArrayListWithExpectedSize(2);
        mPendingReads = Lists.newArrayListWithExpectedSize(2);

        mUnpaired = !read.getReadPairedFlag();

        processRead(read);

        mCachedReads = Lists.newArrayListWithExpectedSize(2);
    }

    public String readId() { return mReadId; }
    public List<ReadInfo> receivedReads() { return mReceivedReads; }
    public List<ReadInfo> pendingReads() { return mPendingReads; }

    public List<SAMRecord> cachedReads() { return mCachedReads; }
    public void addCachedRead(final SAMRecord read) { mCachedReads.add(read); }
    public void clearCachedReads() { mCachedReads.clear(); }

    public boolean processRead(final SAMRecord read)
    {
        // returns true if the read has not been seen before
        ReadInfo readInfo = ReadInfo.fromRead(read);

        if(hasMatch(readInfo, mReceivedReads))
            return false;

        mReceivedReads.add(readInfo);

        for(int i = 0; i < mPendingReads.size(); ++i)
        {
            if(mPendingReads.get(i).matches(readInfo))
            {
                mPendingReads.remove(i);
                break;
            }
        }

        if(read.getReadPairedFlag())
        {
            ReadInfo mateInfo = ReadInfo.fromReadMate(read);

            if(!hasMatch(mateInfo, mReceivedReads) && !hasMatch(mateInfo, mPendingReads))
            {
                mPendingReads.add(mateInfo);
            }
        }

        List<SupplementaryReadData> suppDataList = SupplementaryReadData.extractAlignments(read);

        if(suppDataList != null)
        {
            boolean readIsPrimary = readInfo.IsPrimary;

            for(SupplementaryReadData suppData : suppDataList)
            {
                ReadInfo suppReadInfo = ReadInfo.fromSupplementaryData(suppData, !readIsPrimary, readInfo.FirstInPair);

                if(!hasMatch(suppReadInfo, mReceivedReads) && !hasMatch(suppReadInfo, mPendingReads))
                {
                    mPendingReads.add(suppReadInfo);
                }

                if(!readIsPrimary) // any further supplementaries are about other supplementaries, and so do not need to be recorded
                    break;
            }
        }

        return true;
    }

    private static boolean hasMatch(final ReadInfo readInfo, final List<ReadInfo> existingReadInfos)
    {
        return existingReadInfos.stream().anyMatch(x -> x.matches(readInfo));
    }

    public boolean hasPrimaries()
    {
        if(mReceivedReads.stream().noneMatch(x -> x.IsPrimary && x.FirstInPair))
            return false;

        if(mUnpaired)
            return true;

        return mReceivedReads.stream().anyMatch(x -> x.IsPrimary && !x.FirstInPair);
    }

    public boolean isComplete()
    {
        // only valid if all primaries have been received
        if(!hasPrimaries())
            return false;

        return mPendingReads.isEmpty();
    }

    public List<BasePosition> extractPendingRegions()
    {
        if(mPendingReads.isEmpty())
            return Collections.emptyList();

        List<BasePosition> regions = Lists.newArrayListWithExpectedSize(mPendingReads.size());

        for(ReadInfo readInfo : mPendingReads)
        {
            regions.add(new BasePosition(readInfo.Contig, readInfo.AlignmentStart));
        }

        return regions;
    }

    public String toString()
    {
        return format("id(%s) received(%d) pending(%d) cached(%d)",
                mReadId, mReceivedReads.size(), mPendingReads.size(), mCachedReads.size());
    }
}
