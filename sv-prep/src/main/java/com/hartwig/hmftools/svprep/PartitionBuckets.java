package com.hartwig.hmftools.svprep;

import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;

import java.util.function.Consumer;

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
}
