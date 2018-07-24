package com.hartwig.hmftools.data_analyser.types;

import static java.lang.Math.floor;
import static java.lang.Math.log10;
import static java.lang.Math.pow;
import static java.lang.Math.round;

import java.util.List;

import com.google.common.collect.Lists;

public class BucketGroup implements Comparable<BucketGroup> {

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

    public int getScore() { return mBucketIds.size() * mSampleIds.size(); }

    public int compareTo(final BucketGroup other)
    {
        // for descending order
        return other.getScore() - getScore();
    }

    public final List<Integer> getSampleIds() { return mSampleIds; }
    public final List<Integer> getBucketIds() { return mBucketIds; }
    public final List<BucketPair> getBucketPairs() { return mBucketPairs; }

    public void addBucketPair(final BucketPair bucketPair)
    {
        if(!mBucketPairs.contains(bucketPair))
            mBucketPairs.add(bucketPair);

        // maintain a list of unique sample and bucket IDs
        if(!mBucketIds.contains(bucketPair.getBucketA()))
            mBucketIds.add(bucketPair.getBucketA());

        if(!mBucketIds.contains(bucketPair.getBucketB()))
            mBucketIds.add(bucketPair.getBucketB());
    }

    public boolean hasSample(int sampleIndex)
    {
        return mSampleIds.contains(sampleIndex);
    }

    public void addSamples(List<Integer> sampleIds)
    {
        for(Integer sample : sampleIds)
        {
            addSample(sample);
        }
    }

    public void addSample(int sampleId)
    {
        if(mSampleIds.contains(sampleId))
            return;

        mSampleIds.add(sampleId);
    }

    public boolean hasBucket(int bucketIndex)
    {
        return mBucketIds.contains(bucketIndex);
    }

    public void addBuckets(List<Integer> bucketIds)
    {
        for(Integer bucket : bucketIds)
        {
            addBucket(bucket);
        }
    }

    public void addBucket(int bucketId)
    {
        if(mBucketIds.contains(bucketId))
            return;

        mBucketIds.add(bucketId);
    }

    public static List<Integer> getMatchingBucketList(final List<Integer> bl1, final List<Integer> bl2)
    {
        List<Integer> matchedList = Lists.newArrayList();

        for(Integer bucket : bl1)
        {
            if(bl2.contains(bucket))
                matchedList.add(bucket);
        }

        return matchedList;
    }

    public static boolean hasMatchingBucketList(final List<Integer> bl1, final List<Integer> bl2, double matchPercent)
    {
        int matched = 0;

        for(Integer bucket : bl1)
        {
            if(bl2.contains(bucket))
                ++matched;
        }

        return matched / (double)bl1.size() >= matchPercent;

    }

}
