package com.hartwig.hmftools.data_analyser.types;

import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getSortedVectorIndices;
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
    List<Integer> mSampleIds;

    double[] mCombinedCounts;
    double[] mBucketCountRatios;
    double[] mBucketMeanRatios;
    double[] mBucketRatiosHigh;
    double[] mBucketRatiosLow;
    double mTotalCount;

    private String mCancerType;
    private String mEffects;

    private static final Logger LOGGER = LogManager.getLogger(BucketFamily.class);

    public BucketFamily(int id)
    {
        mId = id;
        mBucketIds = Lists.newArrayList();
        mSampleIds = Lists.newArrayList();
        mBucketGroups = Lists.newArrayList();

        mCancerType = "";
        mEffects = "";
    }

    public int getId() { return mId; }

    public final List<Integer> getBucketIds() { return mBucketIds; }
    public final List<Integer> getSampleIds() { return mSampleIds; }
    public final double[] getBucketCounts() { return mCombinedCounts; }
    public final double[] getBucketRatios() { return mBucketMeanRatios; }
    public final double[] getBucketHighRatios() { return mBucketRatiosHigh; }
    public final double[] getBucketLowRatios() { return mBucketRatiosLow; }
    public double getTotalCount() { return mTotalCount; }

    public List<BucketGroup> getBucketGroups() { return mBucketGroups; }
    public final String getCancerType() { return mCancerType; }
    public final String getEffects() { return mEffects; }

    public void addBucketGroup(BucketGroup group)
    {
        mBucketGroups.add(group);

        for(Integer bucket : group.getBucketIds())
        {
            if(!mBucketIds.contains(bucket))
                mBucketIds.add(bucket);
        }

        for(Integer sample : group.getSampleIds())
        {
            if(!mSampleIds.contains(sample))
                mSampleIds.add(sample);
        }

        // merge cancer type and effects
        if(!group.getCancerType().isEmpty())
        {
            String[] newTypes = group.getCancerType().split(";");

            for(final String type : newTypes)
            {
                if (mCancerType.isEmpty())
                    mCancerType = type;
                else if (!mCancerType.contains(type))
                    mCancerType += ";" + type;
            }
        }

        if(!group.getEffects().isEmpty())
        {
            String[] newTypes = group.getEffects().split(";");

            for(final String type : newTypes)
            {
                if (mEffects.isEmpty())
                    mEffects = type;
                else if (!mEffects.contains(type))
                    mEffects += ";" + type;
            }
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

    public void calcAll(final NmfMatrix sampleCounts)
    {
        if(mBucketGroups.isEmpty())
            return;

        final double[][] scData = sampleCounts.getData();

        // initialise data arrays
        int bucketCount = mBucketGroups.get(0).getBucketCounts().length;

        mCombinedCounts = new double[bucketCount];
        mBucketCountRatios = new double[bucketCount];
        mBucketMeanRatios = new double[bucketCount];
        mBucketRatiosHigh = new double[bucketCount];
        mBucketRatiosLow = new double[bucketCount];

        int[] bgCount = new int[bucketCount];
        int[] sampleCount = new int[bucketCount];

        // work out sample bucket ratios just for the buckets in this family
        List<Double> sampleTotals = Lists.newArrayList();

        for(Integer sampleId : mSampleIds)
        {
            double sampleBucketTotal = 0;

            for (Integer bucket : mBucketIds)
            {
                sampleBucketTotal += sampleCounts.get(bucket, sampleId);
            }

            sampleTotals.add(sampleBucketTotal);
        }

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

            /*
            for(int samIndex = 0; samIndex < mSampleIds.size(); ++samIndex)
            {
                Integer sampleId = mSampleIds.get(samIndex);

                double sbCount = scData[bucket][sampleId];
                double sbRatio = sbCount / sampleTotals.get(samIndex);

                mCombinedCounts[bucket] += sbCount;

                // put the counts for each sample into ascending order so percentiles can be calculated
                int i = 0;
                for (; i < ascRatios.size(); ++i)
                {
                    if (sbRatio < ascRatios.get(i))
                        break;
                }

                ascRatios.add(i, sbRatio);
                ascCounts.add(i, sbCount);
            }
            */

            double bucketTotal = mCombinedCounts[bucket];

            // set the high and low ratios as the 90th and 10th percentiles from the counts
            double lowCountPercentile = 0.25 * bucketTotal;
            double highCountPercentile = 0.75 * bucketTotal;
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

        // log in order of highest ML down
        List<Integer> bucketIndexByRatio = getSortedVectorIndices(mBucketMeanRatios, false);

        for(Integer bucket : bucketIndexByRatio)
        {
            if(mBucketMeanRatios[bucket] == 0)
                break;

            // Integer bucket = mCombinedCounts.get(bucketIndex);

            mBucketCountRatios[bucket] = mCombinedCounts[bucket] / mTotalCount;

            LOGGER.debug(String.format("family(%d) bucket(%d) count(%.0f) ratio(mean=%.3f range=%.3f -> %.3f)",
                    mId, bucket, mCombinedCounts[bucket],
                    mBucketMeanRatios[bucket], mBucketRatiosLow[bucket], mBucketRatiosHigh[bucket]));
        }
    }

}
