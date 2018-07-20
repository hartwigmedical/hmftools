package com.hartwig.hmftools.data_analyser.types;

import static java.lang.Math.floor;
import static java.lang.Math.log10;
import static java.lang.Math.pow;
import static java.lang.Math.round;

import java.util.List;

import com.google.common.collect.Lists;

public class BucketGroup {

    // keyed by a bucket pairing
    int mId;

    List<Integer> mSampleIds;
    List<Integer> mBucketIds;
    List<BucketPair> mBucketPairs;

    public BucketGroup(int id, List<Integer> sampleIds)
    {
        mId = id;

        mSampleIds = sampleIds;
        mBucketIds = Lists.newArrayList();
        mBucketPairs = Lists.newArrayList();
    }

    public int getId() { return mId; }

    public final List<Integer> getSampleIds() { return mSampleIds; }
    public final List<Integer> getBucketIds() { return mBucketIds; }
    public final List<BucketPair> getBucketPairs() { return mBucketPairs; }

    public void addBucketPair(final BucketPair bucketPair)
    {
        mBucketPairs.add(bucketPair);

        // maintain a list of unique sample and bucket IDs
        if(!mBucketIds.contains(bucketPair.getBucketA()))
            mBucketIds.add(bucketPair.getBucketA());

        if(!mBucketIds.contains(bucketPair.getBucketB()))
            mBucketIds.add(bucketPair.getBucketB());
    }

    public void addSamples(List<Integer> sampleIds)
    {
        for(Integer newSample : sampleIds)
        {
            if(!hasSample(newSample))
                mSampleIds.add(newSample);
        }
    }

    public boolean hasSample(int sampleIndex)
    {
        return mSampleIds.contains(sampleIndex);
    }

    public boolean hasBucket(int bucketIndex)
    {
        return mBucketIds.contains(bucketIndex);
    }
}
