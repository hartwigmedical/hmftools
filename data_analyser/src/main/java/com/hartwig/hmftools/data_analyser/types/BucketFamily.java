package com.hartwig.hmftools.data_analyser.types;

import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sumVector;

import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class BucketFamily
{
    private int mId;

    private List<BucketGroup> mBucketGroups;
    List<Integer> mBucketIds;

    double[] mCombinedCounts;
    double[] mBucketCountRatios;
    double[] mBucketMeanRatios;
    double[] mBucketRatiosHigh;
    double[] mBucketRatiosLow;
    double mTotalCount;

    private static final Logger LOGGER = LogManager.getLogger(BucketFamily.class);

    public BucketFamily(int id)
    {
        mId = id;
        mBucketIds = Lists.newArrayList();
        mBucketGroups = Lists.newArrayList();
    }

    public int getId() { return mId; }

    public final List<Integer> getBucketIds() { return mBucketIds; }
    public final double[] getBucketCounts() { return mCombinedCounts; }
    public final double[] getBucketRatios() { return mBucketMeanRatios; }
    public final double[] getBucketHighRatios() { return mBucketRatiosHigh; }
    public final double[] getBucketLowRatios() { return mBucketRatiosLow; }

    public List<BucketGroup> getBucketGroups() { return mBucketGroups; }

    public void addBucketGroup(BucketGroup group)
    {
        mBucketGroups.add(group);

        for(Integer bucket : group.getBucketIds())
        {
            if(!mBucketIds.contains(bucket))
                mBucketIds.add(bucket);
        }
    }

    public boolean hasBucketGroup(final BucketGroup group)
    {
        for(final BucketGroup bucketGroup : mBucketGroups)
        {
            if(group == bucketGroup)
                return true;
        }

        return false;
    }

    public void calcAll()
    {
        if(mBucketGroups.isEmpty())
            return;

        // initialise data arrays
        int bucketCount = mBucketGroups.get(0).getBucketCounts().length;

        mCombinedCounts = new double[bucketCount];
        mBucketCountRatios = new double[bucketCount];
        mBucketMeanRatios = new double[bucketCount];
        mBucketRatiosHigh = new double[bucketCount];
        mBucketRatiosLow = new double[bucketCount];

        int[] bgCount = new int[bucketCount];
        int[] sampleCount = new int[bucketCount];

        for (Integer bucket : mBucketIds)
        {
            List<Double> ascRatios = Lists.newArrayList();
            List<Double> ascCounts = Lists.newArrayList();

            for(final BucketGroup bucketGroup : mBucketGroups)
            {
                if(!bucketGroup.hasBucket(bucket))
                    continue;

                final double[] bucketRatios = bucketGroup.getBucketRatios();
                final double[] bucketCounts = bucketGroup.getBucketCounts();

                mCombinedCounts[bucket] += bucketCounts[bucket];

                ++bgCount[bucket];
                sampleCount[bucket] += bucketGroup.getSampleIds().size();

                // put the ratios for each bucket into ascending order so percentiles can be calculated
                int i = 0;
                for (; i < ascRatios.size(); ++i)
                {
                    if (bucketRatios[bucket] < ascRatios.get(i))
                        break;
                }

                ascRatios.add(i, bucketRatios[bucket]);
                ascCounts.add(i, bucketCounts[bucket]);
            }

            double bucketTotal = mCombinedCounts[bucket];

            // set the high and low ratios as the 90th and 10th percentiles from the counts
            double lowCountPercentile = 0.1 * bucketTotal;
            double highCountPercentile = 0.9 * bucketTotal;
            double workingCount = 0;
            double workingRatioTotal = 0;
            double workingCountLow = 0;
            double workingRatioTotalLow = 0;
            double workingCountHigh = 0;
            double workingRatioTotalHigh = 0;

            for(int i = 0; i < ascRatios.size(); ++i)
            {
                double ratio = ascRatios.get(i);
                double count = ascCounts.get(i);

                workingCount += count;
                workingRatioTotal += count * ratio;

                if(workingCountLow + count >= lowCountPercentile)
                {
                    if(workingCountLow < lowCountPercentile)
                    {
                        workingRatioTotalLow += (lowCountPercentile - workingCountLow) * ratio;
                        workingCountLow = lowCountPercentile;
                    }
                }
                else
                {
                    workingCountLow += count;
                    workingRatioTotalLow += count * ratio;
                }

                if(workingCountHigh == 0)
                {
                    // haven't yet high the high-percentile mark
                    if(workingCount > highCountPercentile)
                    {
                        workingCountHigh = workingCount - highCountPercentile;
                        workingRatioTotalHigh = workingCountHigh * ratio;
                    }
                }
                else
                {
                    workingCountHigh += count;
                    workingRatioTotalHigh += count * ratio;
                }
            }

            mBucketRatiosLow[bucket] = workingRatioTotalLow / workingCountLow;
            mBucketRatiosHigh[bucket] = workingRatioTotalHigh / workingCountHigh;
            mBucketMeanRatios[bucket] = workingRatioTotal / workingCount;
        }

        mTotalCount = sumVector(mCombinedCounts);

        // form bucket ratios the old way
        for (Integer bucket : mBucketIds)
        {
            mBucketCountRatios[bucket] = mCombinedCounts[bucket] / mTotalCount;

            LOGGER.debug(String.format("bf(%d) bucket(%d groups=%d samples=%d) count(%.0f) ratio(mean=%.3f range=%.3f -> %.3f)",
                    mId, bucket, bgCount[bucket], sampleCount[bucket], mCombinedCounts[bucket],
                    mBucketMeanRatios[bucket], mBucketRatiosLow[bucket], mBucketRatiosHigh[bucket]));
        }
    }

}
