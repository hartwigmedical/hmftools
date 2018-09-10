package com.hartwig.hmftools.data_analyser.types;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.data_analyser.calcs.BucketAnalyser.MAX_NOISE_ALLOC_PERC;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.capValue;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.copyVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.greaterThan;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.initVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.lessThan;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sumVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.vectorMultiply;
import static com.hartwig.hmftools.data_analyser.types.BucketGroup.ratioRange;

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
    private double[] mAllocNoiseCounts;
    private double[] mAllocBucketCounts;
    private double[] mUnallocBucketCounts;
    private double[] mPartialUnallocBucketCounts; // purely for the discovery phase, holding some allocated counts back
    private double mAllocPercent;
    private double mPreviousAllocPerc; // to allow tracking of changes
    private double mVarTotal;
    private double mElevatedTotal;
    private double mAllocTotal;
    private double mUnallocTotal;
    private double mNoiseAllocTotal;
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
        mPreviousAllocPerc = 0;
        mAllocTotal = 0;
        mUnallocTotal = 0;
        mNoiseAllocTotal = 0;
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

    public final BucketGroup getBackgroundGroup() { return mBackgroundGroup; }

    public void setBackgroundGroup(final BucketGroup group)
    {
        mBackgroundGroup = group;
    }

    public void addElevBucketGroup(final BucketGroup group, double allocPerc)
    {
        mElevBucketGroups.add(group);
        mGroupAllocPercents.add(allocPerc);
    }

    public final double[] getElevatedBucketCounts() { return mElevBucketCounts; }
    public final double[] getCountRanges() { return mCountRanges; }
    public final double[] getAllocNoiseCounts() { return mAllocNoiseCounts; }
    public double getTotalCount() { return mVarTotal; }
    public double getElevatedCount() { return mElevatedTotal; }
    public double getAllocatedCount() { return mAllocTotal; }
    public double getUnallocatedCount() { return mUnallocTotal; }
    public double getAllocNoise() { return mNoiseAllocTotal; }
    public double getNoiseTotal() { return mNoiseTotal; }

    public final double[] getUnallocBucketCounts() { return mUnallocBucketCounts; }
    public final double[] getPartialUnallocBucketCounts() { return mPartialUnallocBucketCounts; }
    public double getAllocPercent() { return mAllocPercent; }
    public double lastAllocPercChange() { return mAllocPercent - mPreviousAllocPerc; }
    public double getUnallocPercent() { return 1 - mAllocPercent; }
    public double getNoisePerc() { return mNoiseAllocTotal/mNoiseTotal; }
    public double getNoiseOfTotal() { return mNoiseAllocTotal/mElevatedTotal; }

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
        mAllocNoiseCounts = new double[counts.length];
        mUnallocBucketCounts = new double[counts.length];
        mPartialUnallocBucketCounts = new double[counts.length];
        mCountRanges = new double[counts.length];
        copyVector(counts, mUnallocBucketCounts);
        copyVector(counts, mPartialUnallocBucketCounts);
        copyVector(counts, mElevBucketCounts);
        mElevatedTotal = sumVector(counts);
        mUnallocTotal = mElevatedTotal;

        copyVector(ranges, mCountRanges);
        mNoiseTotal = sumVector(mCountRanges);
    }

    public void clearAllocations()
    {
        mElevBucketGroups.clear();
        mUnallocBuckets.clear();
        mUnallocBuckets.addAll(mElevatedBuckets);

        copyVector(mElevBucketCounts, mUnallocBucketCounts);
        copyVector(mElevBucketCounts, mPartialUnallocBucketCounts);
        initVector(mAllocBucketCounts, 0);
        initVector(mAllocNoiseCounts, 0);
        mUnallocTotal = mElevatedTotal;
        mGroupAllocPercents.clear();
        mAllocPercent = 0;
        mPreviousAllocPerc = 0;
        mAllocTotal = 0;
        mNoiseAllocTotal = 0;
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

    public double[] getPotentialUnallocCounts(final double[] bucketRatios, final List<Integer> requiredBuckets, final double[] ratioRanges)
    {
        return getPotentialAllocation(bucketRatios, true, requiredBuckets, ratioRanges);
    }

    public double[] getPotentialElevCounts(final double[] bucketRatios, final List<Integer> requiredBuckets, final double[] ratioRanges)
    {
        return getPotentialAllocation(bucketRatios, false, requiredBuckets, ratioRanges);
    }

    private double[] getPotentialAllocation(final double[] bucketRatios, boolean useUnallocated, final List<Integer> requiredBuckets, final double[] ratioRanges)
    {
        // must cap at the actual sample counts
        // allow to go as high as the elevated probability range
        // if a single bucket required by the ratios has zero unallocated, then all are zeroed
        final double[] sampleCounts = useUnallocated ? mUnallocBucketCounts : mElevBucketCounts;
        double countsTotal = useUnallocated ? mUnallocTotal : mElevatedTotal;

        // first extract the remaining unallocated counts per bucket
        double minAlloc = 0;

        for(int bIndex = 0; bIndex < requiredBuckets.size(); ++bIndex)
        {
            Integer bucket = requiredBuckets.get(bIndex);
            double noiseCount = useUnallocated ? max(mCountRanges[bucket] - mAllocNoiseCounts[bucket], 0) : mCountRanges[bucket];

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
            double potentialAlloc = unallocCount + noiseCount;

            // use the low range ratio to give max possible allocation to this bucket
            double ratio = bucketRatios[bucket] + ratioRange(ratioRanges, bucket, true);
            double alloc = potentialAlloc / ratio;

            if(minAlloc == 0 || alloc < minAlloc)
                minAlloc = alloc;
        }

        // cap by the actual unallocated before working out allocation
        minAlloc = capValue(minAlloc, 0, countsTotal);

        double[] bestAllocCounts = new double[bucketRatios.length];

        if(minAlloc == 0)
            return bestAllocCounts;

        double noiseTotal = useUnallocated ? mNoiseAllocTotal : 0;

        for(int bIndex = 0; bIndex < requiredBuckets.size(); ++bIndex)
        {
            Integer bucket = requiredBuckets.get(bIndex);

            double ratio = bucketRatios[bucket] + ratioRange(ratioRanges, bucket, false);
            double potentialAlloc = minAlloc * ratio;

            if(potentialAlloc == 0)
                continue;

            double noiseCount = useUnallocated ? max(mCountRanges[bucket] - mAllocNoiseCounts[bucket], 0) : mCountRanges[bucket];

            if(noiseCount + noiseTotal > mVarTotal * MAX_NOISE_ALLOC_PERC)
            {
                noiseCount = max(mVarTotal * MAX_NOISE_ALLOC_PERC - noiseTotal, 0);
            }

            double allocCount = potentialAlloc;

            if(greaterThan(allocCount, sampleCounts[bucket] + noiseCount))
                allocCount = sampleCounts[bucket] + noiseCount; // shouldn't happen since should have determined above

            bestAllocCounts[bucket] = allocCount;

            if(noiseCount > 0 && sampleCounts[bucket] < allocCount)
            {
                noiseTotal += allocCount - sampleCounts[bucket];
            }
        }

        return bestAllocCounts;
    }

    public double allocateBucketCounts(double[] counts, double reqAllocationPercent)
    {
        double allocatedCount = 0;
        double allocatedActualCount = 0;
        double noiseTotal = mNoiseAllocTotal;
        double reductionFactor = 1;

        // do a preliminary check that this allocation will actually achieve the required percentage count
        for(int i = 0; i < counts.length; ++i)
        {
            if(counts[i] == 0)
                continue;

            double unallocNoise = max(mCountRanges[i] - mAllocNoiseCounts[i], 0);

            if(unallocNoise + noiseTotal > mVarTotal * MAX_NOISE_ALLOC_PERC)
            {
                unallocNoise = max(mVarTotal * MAX_NOISE_ALLOC_PERC - noiseTotal, 0);
            }

            double allocCount = min(mUnallocBucketCounts[i] + unallocNoise, counts[i]);

            reductionFactor = min(reductionFactor, allocCount / counts[i]);

            // if even a single bucket restricts a non-zero count being allocated, fail the whole allocation
            if(allocCount == 0)
                return 0;

            allocatedCount += allocCount;

            if(unallocNoise > 0 && mUnallocBucketCounts[i] < allocCount)
            {
                noiseTotal += allocCount - mUnallocBucketCounts[i];
            }

            if(mUnallocBucketCounts[i] < counts[i])
            {
                allocatedActualCount += mUnallocBucketCounts[i];
                // allocatedNoise += min(counts[i] - mUnallocBucketCounts[i], unallocNoise);
            }
            else
            {
                allocatedActualCount += counts[i];
            }
        }

        if(lessThan(reductionFactor, 1))
        {
            // assumption is that the counts are in proportion (eg if generated off
            allocatedActualCount *= reductionFactor;
            allocatedCount *= reductionFactor;
            vectorMultiply(counts, reductionFactor);
        }

        // now only factors in the allocation to remaining actual counts
        if(reqAllocationPercent > 0 && allocatedCount < reqAllocationPercent * mElevatedTotal)
            return allocatedCount;

        // allow allocation for the caller to go up to the unallocated counts plus the permitted noise range, which
        // will only be allocated once (when count > unallocated)
        // modify the caller's count values if limited in any way
        // internal allocation counts and total are limited to actuals
        noiseTotal = mNoiseAllocTotal;

        for(int i = 0; i < counts.length; ++i)
        {
            if(counts[i] == 0)
                continue;

            if (mUnallocBucketCounts[i] < counts[i])
            {
                double unallocNoise = max(mCountRanges[i] - mAllocNoiseCounts[i], 0);

                if(unallocNoise + noiseTotal > mVarTotal * MAX_NOISE_ALLOC_PERC)
                {
                    unallocNoise = max(mVarTotal * MAX_NOISE_ALLOC_PERC - noiseTotal, 0);
                }

                if(mUnallocBucketCounts[i] == 0 && unallocNoise == 0)
                {
                    counts[i] = 0;
                    continue;
                }

                if(unallocNoise > 0)
                {
                    // if unalloc = 10, noise = 15, count = 20, then just allocate 10 of noise
                    // if unalloc = 10, noise = 5, count = 20, then allocate all 5 of noise
                    double noiseAlloc = min(counts[i] - mUnallocBucketCounts[i], unallocNoise);

                    mNoiseAllocTotal += noiseAlloc; // eg if unalloc = 10, noise = 5, count = 20, then allocate all 5 of noise
                    mAllocNoiseCounts[i] = noiseAlloc;
                    noiseTotal += noiseAlloc;
                }

                if(mUnallocBucketCounts[i] == 0)
                {
                    // just check if any remaining noise needs to be allocated
                    counts[i] = unallocNoise;
                    noiseTotal += unallocNoise;
                    continue;
                }
                else
                {
                    mAllocBucketCounts[i] += mUnallocBucketCounts[i];

                    counts[i] = min(mUnallocBucketCounts[i] + unallocNoise, counts[i]);
                    mUnallocBucketCounts[i] = 0;
                    removeUnallocBucket(i);
                }
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

        mPreviousAllocPerc = mAllocPercent;
        mAllocPercent = capValue(mAllocTotal/mElevatedTotal, 0, 1);

        return allocatedCount;
    }

    private void removeUnallocBucket(Integer bucketId)
    {
        mUnallocBuckets.remove(bucketId);
    }

}
