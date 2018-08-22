package com.hartwig.hmftools.data_analyser.types;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.capValue;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.copyVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.initVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sumVector;

import java.util.List;

import com.google.common.collect.Lists;

public class SampleData
{
    final public int Id;

    private String mSampleName;
    private boolean mExcluded;

    private double[] mBucketCounts;
    private double[] mBackgroundCounts; // those assigned to the background signature
    private double[] mElevBucketCounts;
    private double[] mCountRanges;
    private double[] mAllocBucketCounts;
    private double[] mUnallocBucketCounts;
    private double[] mPartialUnallocBucketCounts; // purely for the discovery phase, holding some allocated counts back
    private double mAllocPercent;
    private double mVarTotal;
    private double mElevatedTotal;
    private double mAllocTotal;
    private double mUnallocTotal;
    private double mNoiseAlloc;
    private double mNoiseTotal;

    private String mCancerType;

    private final List<Integer> mElevatedBuckets;
    private final List<Integer> mUnallocBuckets;
    private final List<String> mCategoryData;
    private List<BucketGroup> mElevBucketGroups;
    private BucketGroup mBackgroundGroup;
    private final List<Double> mGroupAllocPercents;

    public static double PARTIAL_ALLOC_PERCENT = 0.2;

    public SampleData(int id)
    {
        Id = id;
        mSampleName = "";
        mExcluded = false;
        mCancerType = "";
        mCategoryData = Lists.newArrayList();
        mElevBucketGroups = Lists.newArrayList();
        mElevatedBuckets = Lists.newArrayList();
        mUnallocBuckets = Lists.newArrayList();
        mGroupAllocPercents = Lists.newArrayList();
        mBackgroundGroup = null;
        mAllocPercent = 0;
        mAllocTotal = 0;
        mUnallocTotal = 0;
        mNoiseAlloc = 0;
        mNoiseTotal = 0;
        mVarTotal = 0;
    }

    public final String getSampleName() { return mSampleName; }
    public void setSampleName(final String name) { mSampleName = name; }

    public final String getCancerType() { return mCancerType; }
    public void setCancerType(final String type) { mCancerType = type; }

    public final List<String> getCategoryData() { return mCategoryData; }
    public void setCategoryData(final List<String> data) { mCategoryData.addAll(data); }

    public final boolean isExcluded() { return mExcluded; }
    public void setExcluded(boolean toggle) { mExcluded = toggle; }

    public final void setElevatedBuckets(List<Integer> buckets)
    {
        mElevatedBuckets.addAll(buckets);
        mUnallocBuckets.addAll(buckets);
    }

    public final List<Integer> getElevatedBuckets() { return mElevatedBuckets; }
    public final List<Integer> getUnallocBuckets() { return mUnallocBuckets; }

    public final List<BucketGroup> getElevBucketGroups() { return mElevBucketGroups; }
    public final List<Double> getGroupAllocPercents() { return mGroupAllocPercents; }
    public void clearElevBucketGroups()
    {
        mElevBucketGroups.clear();
        mGroupAllocPercents.clear();
    }

    public final BucketGroup getBackgroundGroup() { return mBackgroundGroup; }

    public void setBackgroundGroup(final BucketGroup group)
    {
        mBackgroundGroup = group;
    }

    public void addElevBucketGroup(final BucketGroup group)
    {
        mElevBucketGroups.add(group);
    }

    public void addElevBucketGroup(final BucketGroup group, double allocPerc)
    {
        mElevBucketGroups.add(group);
        mGroupAllocPercents.add(allocPerc);
    }

    public final double[] getBucketCounts() { return mBucketCounts; }
    public final double[] getElevatedBucketCounts() { return mElevBucketCounts; }
    public final double[] getBackgroundCounts() { return mBackgroundCounts; }
    public final double[] getCountRanges() { return mCountRanges; }
    public double getTotalCount() { return mVarTotal; }
    public double getElevatedCount() { return mElevatedTotal; }
    public double getAllocatedCount() { return mAllocTotal; }
    public double getUnallocatedCount() { return mUnallocTotal; }
    public double getAllocNoise() { return mNoiseAlloc; }
    public double getNoiseTotal() { return mNoiseTotal; }

    public final double[] getAllocBucketCounts() { return mAllocBucketCounts; }
    public final double[] getUnallocBucketCounts() { return mUnallocBucketCounts; }
    public final double[] getPartialUnallocBucketCounts() { return mPartialUnallocBucketCounts; }
    public double getAllocPercent() { return mAllocPercent; }
    public double getUnallocPercent() { return 1 - mAllocPercent; }

    public void setBucketCounts(final double[] counts)
    {
        mBucketCounts = new double[counts.length];
        copyVector(counts, mBucketCounts);
        mVarTotal = sumVector(counts);
    }

    public void setElevatedBucketCounts(final double[] counts, final double[] ranges)
    {
        mElevBucketCounts = new double[counts.length];
        mAllocBucketCounts = new double[counts.length];
        mUnallocBucketCounts = new double[counts.length];
        mPartialUnallocBucketCounts = new double[counts.length];
        mCountRanges = new double[counts.length];
        copyVector(counts, mUnallocBucketCounts);
        copyVector(counts, mPartialUnallocBucketCounts);
        copyVector(counts, mElevBucketCounts);
        copyVector(ranges, mCountRanges);
        mElevatedTotal = sumVector(counts);
        mUnallocTotal = mElevatedTotal;
        mNoiseTotal = sumVector(ranges);
    }

    public void clearAllocations()
    {
        mElevBucketGroups.clear();
        mUnallocBuckets.clear();
        mUnallocBuckets.addAll(mElevatedBuckets);

        copyVector(mElevBucketCounts, mUnallocBucketCounts);
        copyVector(mElevBucketCounts, mPartialUnallocBucketCounts);
        mAllocBucketCounts = new double[mUnallocBucketCounts.length];
        mUnallocTotal = mElevatedTotal;
        mGroupAllocPercents.clear();
        mAllocPercent = 0;
        mAllocTotal = 0;
        mNoiseAlloc = 0;
    }

    public void populateBucketCountSubset(double[] counts, final List<Integer> bucketSubset, boolean usePartials)
    {
        // extract the counts for the specified subset, leaving the rest zeroed
        initVector(counts, 0);

        for(Integer bucketId : bucketSubset)
        {
            if(usePartials)
                counts[bucketId] = max(mPartialUnallocBucketCounts[bucketId], 0);
            else
            counts[bucketId] = max(mUnallocBucketCounts[bucketId], 0);
        }
    }

    public double[] getPotentialUnallocCounts(final double[] bucketRatios, final List<Integer> requiredBuckets)
    {
        return getPotentialAllocation(bucketRatios, mUnallocBucketCounts, mUnallocTotal, requiredBuckets, null);
    }

    public double[] getPotentialElevCounts(final double[] bucketRatios, final List<Integer> requiredBuckets)
    {
        return getPotentialAllocation(bucketRatios, mElevBucketCounts, mElevatedTotal, requiredBuckets, null);
    }

    public double[] getPotentialAllocation(final double[] bucketRatios, final List<Integer> requiredBuckets, final List<Double> bucketRatioRanges)
    {
        return getPotentialAllocation(bucketRatios, mUnallocBucketCounts, mUnallocTotal, requiredBuckets, bucketRatioRanges);
    }

    public double[] getPotentialAllocation(final double[] bucketRatios, final double[] sampleCounts, double countsTotal,
            final List<Integer> requiredBuckets, final List<Double> bucketRatioRanges)
    {
        // must cap at the actual sample counts
        // allow to go as high as the elevated probability range
        // if a single bucket required by the ratios has zero unallocated, then all are zeroed

        // first extract the remaining unallocated counts per bucket
        double[] bestAllocCounts = new double[bucketRatios.length];

        // TEMP: ratio rangs not used
        final List<Double> ratioRanges = null;

        double minAlloc = 0;

        for(int bIndex = 0; bIndex < requiredBuckets.size(); ++bIndex)
        {
            Integer bucket = requiredBuckets.get(bIndex);
            double ratioRange = ratioRanges != null ? ratioRanges.get(bIndex) : 0;

            // if any of the required buckets are already fully allocated, no others can be allocated
            double unallocCount = max(sampleCounts[bucket], 0);

            if(unallocCount == 0)
            {
                if(mBucketCounts[bucket] > 0 && mElevBucketCounts[bucket] > 0)
                {
                    // fully allocated
                    minAlloc = 0;
                    break;
                }

                // otherwise it is valid to use non-negative noise with a zero elevated count
            }

            // allow full allocation of noise, not proportional to allocated counts to make consistent with fitter
            // double unallocPerc = unallocCount / mElevBucketCounts[bucket];
            double potentialAlloc = unallocCount + mCountRanges[bucket]; // unallocPerc * mCountRanges[bucket];

            // use the low range ratio to give max possible allocation to this bucket
            double adjBucketRatio = bucketRatios[bucket] * (1 - ratioRange);

            double alloc = potentialAlloc / adjBucketRatio;

            if(minAlloc == 0 || alloc < minAlloc)
                minAlloc = alloc;
        }

        // cap by the actual unallocated before working out allocation
        minAlloc = capValue(minAlloc, 0, countsTotal);

        for(int bIndex = 0; bIndex < requiredBuckets.size(); ++bIndex)
        {
            Integer bucket = requiredBuckets.get(bIndex);
            double ratioRange = ratioRanges != null ? ratioRanges.get(bIndex) : 0;

            double potentialAlloc = minAlloc * bucketRatios[bucket] * (1 + ratioRange);
            bestAllocCounts[bucket] = min(sampleCounts[bucket] + mCountRanges[bucket], potentialAlloc);
        }

        return bestAllocCounts;
    }

    public double allocateBucketCounts(double[] counts, double reqAllocationPercent)
    {
        double allocatedCount = 0;

        // do a preliminary check that this allocation will actually apply the required percentage count
        for(int i = 0; i < counts.length; ++i)
        {
            allocatedCount += min(mUnallocBucketCounts[i] + mCountRanges[i], counts[i]);
        }

        if(reqAllocationPercent > 0 && allocatedCount < reqAllocationPercent * mElevatedTotal)
            return allocatedCount;

        // allow allocation up to the unallocated counts plus the permitted noise range, which
        // will only be allocated once (when count > unallocated)
        // modify the caller's count values if limited in any way
        for(int i = 0; i < counts.length; ++i)
        {
            if(mUnallocBucketCounts[i] == 0)
            {
                counts[i] = 0;
                continue;
            }

            if (mUnallocBucketCounts[i] < counts[i])
            {
                mAllocBucketCounts[i] += mUnallocBucketCounts[i];

                if(mUnallocBucketCounts[i] + mCountRanges[i] > counts[i]) // eg if unalloc = 10, noise = 15, count = 20, then just allocate 10 of noise
                    mNoiseAlloc += counts[i] - mUnallocBucketCounts[i];
                else
                    mNoiseAlloc += mCountRanges[i]; // eg if unalloc = 10, noise = 5, count = 20, then allocate all 5 of noise

                counts[i] = min(mUnallocBucketCounts[i] + mCountRanges[i], counts[i]);
                mUnallocBucketCounts[i] = 0;
                removeUnallocBucket(i);
            }
            else
            {
                mUnallocBucketCounts[i] -= counts[i];
                mAllocBucketCounts[i] += counts[i];
            }

            // maintain an extra unallocated count for discovery purposes - eg 10% higher at all times
            mPartialUnallocBucketCounts[i] = min(mUnallocBucketCounts[i] + PARTIAL_ALLOC_PERCENT * mElevBucketCounts[i], mElevBucketCounts[i]);
        }

        mAllocTotal = capValue(sumVector(mAllocBucketCounts), 0, mElevatedTotal);
        mUnallocTotal = mElevatedTotal - mAllocTotal;
        mAllocPercent = capValue(mAllocTotal/mElevatedTotal, 0, 1);
        return allocatedCount;
    }

    private void removeUnallocBucket(Integer bucketId)
    {
        mUnallocBuckets.remove(bucketId);
    }

    public double calcBucketRange(int bucketId)
    {
        if(mElevBucketCounts[bucketId] == 0)
            return 0;

        return mCountRanges[bucketId] / mElevBucketCounts[bucketId];
    }
}
