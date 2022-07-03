package com.hartwig.hmftools.svprep;

import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;

import java.util.Arrays;
import java.util.List;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import htsjdk.samtools.SAMRecord;

public class PartitionBuckets
{
    private final ChrBaseRegion mRegion;

    private final int mBucketSize;
    private final SvBucket[] mBuckets;
    private int mProcessedBucketIndex;

    public PartitionBuckets(final ChrBaseRegion region, int partitionSize, int bucketSize)
    {
        mRegion = region;

        mBucketSize = bucketSize;
        int bucketCount = partitionSize / mBucketSize + 1;
        mBuckets = new SvBucket[bucketCount];
        mProcessedBucketIndex = -1;
    }

    public int getBucketCount()
    {
        return (int)Arrays.stream(mBuckets).filter(x -> x != null).count();
    }

    public SvBucket findBucket(int position)
    {
        int positionOffset = position - mRegion.start();

        int bucketId = positionOffset / mBucketSize;

        if(bucketId < 0 || bucketId >= mBuckets.length)
        {
            SV_LOGGER.error("partition({}) invalid bucket index({}) from position({} offset={}) vs max({})",
                    mRegion, bucketId, position, positionOffset, mBuckets.length);
            return null;
        }

        if(mBuckets[bucketId] == null)
        {
            int bucketPosStart = mRegion.start() + bucketId * mBucketSize;
            int bucketPosEnd = bucketPosStart + mBucketSize - 1;
            ChrBaseRegion bucketRegion = new ChrBaseRegion(mRegion.Chromosome, bucketPosStart, bucketPosEnd);
            mBuckets[bucketId] = new SvBucket(bucketId, bucketRegion);
        }

        return mBuckets[bucketId];
    }

    public void processBuckets(int currentReadPosition, final Consumer<SvBucket> consumer)
    {
        while(mProcessedBucketIndex < mBuckets.length - 1)
        {
            int nextBucketIndex = mProcessedBucketIndex + 1;
            int bucketStartPos = nextBucketIndex * mBucketSize;

            if(currentReadPosition > 0 && currentReadPosition < bucketStartPos - mBucketSize)
                break;

            mProcessedBucketIndex = nextBucketIndex;
            SvBucket bucket = mBuckets[mProcessedBucketIndex];

            if(bucket == null)
                continue;

            consumer.accept(bucket);
        }
    }

    public void transferToNext(final SvBucket bucket)
    {
        if(bucket.id() >= mBuckets.length - 1) // no transfer between partitions
            return;

        SvBucket nextBucket = null;

        // work in descending order until the junctions fall in the previous bucket
        while(!bucket.junctions().isEmpty())
        {
            int index = bucket.junctions().size() - 1;
            JunctionData junctionData = bucket.junctions().get(index);

            if(junctionData.Position <= bucket.region().end())
                break;

            if(nextBucket == null)
                nextBucket = findBucket(junctionData.Position);

            nextBucket.junctions().add(0, junctionData); // add to start since handling in descending position order
            bucket.junctions().remove(index);
        }

        // take any left over support reads, which will straddle this upper bucket
        if(!bucket.supportingReads().isEmpty())
        {
            if(nextBucket == null)
                nextBucket = findBucket(bucket.region().end() + 1);

            nextBucket.supportingReads().addAll(bucket.supportingReads());
            bucket.supportingReads().clear();
        }
    }

    @VisibleForTesting
    public List<SvBucket> getBuckets()
    {
        return Arrays.stream(mBuckets).filter(x -> x != null).collect(Collectors.toList());
    }
}
