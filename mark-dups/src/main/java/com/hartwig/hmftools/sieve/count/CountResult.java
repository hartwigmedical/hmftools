package com.hartwig.hmftools.sieve.count;

import java.util.List;
import java.util.Map;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMSequenceRecord;

public class CountResult
{
    private final SAMSequenceRecord mSeq;
    private final long mPrimaryReadCount;
    private final List<Map<Integer, Long>> mPrimaryBucketCounts;
    private final long mSuppReadCount;
    private final List<Map<Integer, Long>> mSuppBucketCounts;

    public CountResult(
            final SAMSequenceRecord seq,
            final long primaryReadCount,
            @Nullable final List<Map<Integer, Long>> primaryBucketCounts,
            final long suppReadCount,
            @Nullable final List<Map<Integer, Long>> suppBucketCounts)
    {
        mSeq = seq;
        mPrimaryReadCount = primaryReadCount;
        mPrimaryBucketCounts = primaryBucketCounts;
        mSuppReadCount = suppReadCount;
        mSuppBucketCounts = suppBucketCounts;
    }

    public SAMSequenceRecord getSeq()
    {
        return mSeq;
    }

    public long getPrimaryReadCount()
    {
        return mPrimaryReadCount;
    }

    @Nullable
    public List<Map<Integer, Long>> getPrimaryBucketCounts()
    {
        return mPrimaryBucketCounts;
    }

    public long getSuppReadCount()
    {
        return mSuppReadCount;
    }

    @Nullable
    public List<Map<Integer, Long>> getSuppBucketCounts()
    {
        return mSuppBucketCounts;
    }
}
