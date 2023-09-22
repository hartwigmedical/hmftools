package com.hartwig.hmftools.bamtools.slice;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;

import htsjdk.samtools.SAMRecord;

public class ReadGroup
{
    private final String mId;
    private int mSuppCount;
    private int mNonSuppCount;
    private int mExpectedSuppCount;
    private int mExpectedNonSuppCount;

    private boolean mLocalReadsReceived;

    public ReadGroup(final String readId)
    {
        mId = readId;

        mSuppCount = 0;
        mNonSuppCount = 0;
        mExpectedSuppCount = 0;
        mExpectedNonSuppCount = 1;
        mLocalReadsReceived = false;
    }

    public final String id() { return mId; }
    public boolean localReadsReceived() { return mLocalReadsReceived; }

    public void registerMate() { mExpectedNonSuppCount = 2; }

    public void registerSupplementaryData(final SAMRecord read, final SupplementaryReadData suppData)
    {
        if(!read.getSupplementaryAlignmentFlag())
        {
            ++mExpectedSuppCount;
        }
    }

    public void registerRead(final SAMRecord read)
    {
        if(read.getSupplementaryAlignmentFlag())
        {
            ++mSuppCount;
        }
        else
        {
            ++mNonSuppCount;
        }

        mLocalReadsReceived = (mExpectedNonSuppCount == mNonSuppCount) && (mExpectedSuppCount == mSuppCount);
    }

    public String toString()
    {
        return String.format("id(%s) nonSupp(exp=%d rec=%d) supp(exp=%d rec=%d)",
                mId, mExpectedNonSuppCount, mNonSuppCount, mExpectedSuppCount, mSuppCount);
    }
}
