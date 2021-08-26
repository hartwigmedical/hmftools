package com.hartwig.hmftools.sigs.buckets;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.VectorUtils.copyVector;
import static com.hartwig.hmftools.common.utils.VectorUtils.initVector;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;
import static com.hartwig.hmftools.common.utils.VectorUtils.vectorMultiply;
import static com.hartwig.hmftools.sigs.buckets.BaConfig.MAX_NOISE_ALLOC_PERCENT;
import static com.hartwig.hmftools.common.sigs.DataUtils.capValue;
import static com.hartwig.hmftools.common.sigs.DataUtils.greaterThan;
import static com.hartwig.hmftools.common.sigs.DataUtils.lessThan;
import static com.hartwig.hmftools.sigs.buckets.BucketGroup.ratioRange;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SampleData
{
    final public int Id;

    private String mSampleName;
    private boolean mExcluded;

    private double[] mBucketCounts; // counts by bucket
    private double[] mElevBucketCounts; // counts of elevated buckets only
    private double[] mNoiseCounts; // permitted noise or each bucket

    private double[] mAllocBucketCounts; // allocated counts for elevated buckets (excluding noise)
    private double[] mAllocNoiseCounts; // counts allocated to noise for each bucket
    private double[] mUnallocBucketCounts; // unallocated counts for elevated buckets

    private boolean mUseElevatedForAllocation; // any calculated allocations will by from elevated, not total count
    private double mAllocPercent;
    private double mPreviousAllocPerc; // to allow tracking of changes
    private double mVarTotal;
    private double mBackgroundTotal;
    private double mElevatedTotal;
    private double mAllocTotal;
    private double mUnallocTotal;
    private double mMaxNoiseTotal;
    private double mNoiseAllocTotal;

    private String mCancerType;

    private final List<Integer> mElevatedBuckets;
    private final List<Integer> mUnallocBuckets; // of the elevated buckets
    private final List<String> mCategoryData;
    private List<BucketGroup> mBucketGroups;
    private List<BucketGroup> mElevBucketGroups; // doesn't include any background groups
    private BucketGroup mBackgroundGroup;
    private final List<Double> mGroupAllocPercents; // purely informational

    private static final Logger LOGGER = LogManager.getLogger(SampleData.class);

    public SampleData(int id)
    {
        Id = id;
        mSampleName = "";
        mExcluded = false;
        mCancerType = "";
        mCategoryData = Lists.newArrayList();
        mElevBucketGroups = Lists.newArrayList();
        mBucketGroups = Lists.newArrayList();
        mElevatedBuckets = Lists.newArrayList();
        mUnallocBuckets = Lists.newArrayList();
        mGroupAllocPercents = Lists.newArrayList();
        mBackgroundGroup = null;
        mAllocPercent = 0;
        mUseElevatedForAllocation = true;
        mPreviousAllocPerc = 0;
        mAllocTotal = 0;
        mElevatedTotal = 0;
        mBackgroundTotal = 0;
        mUnallocTotal = 0;
        mNoiseAllocTotal = 0;
        mVarTotal = 0;
        mMaxNoiseTotal = 0;
    }

    public final String name() { return mSampleName; }
    public void setSampleName(final String name) { mSampleName = name; }

    public final String cancerType() { return mCancerType; }
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

    public final List<BucketGroup> getBucketGroups() { return mBucketGroups; }
    public final List<BucketGroup> getElevBucketGroups() { return mElevBucketGroups; }

    public final BucketGroup getBackgroundGroup() { return mBackgroundGroup; }

    public final List<Double> getGroupAllocPercents() { return mGroupAllocPercents; }

    public void addBucketGroup(final BucketGroup group, double allocPerc)
    {
        if(mBucketGroups.contains(group))
            return;

        mBucketGroups.add(group);

        if(group.isBackground())
            mBackgroundGroup = group;
        else
            mElevBucketGroups.add(group);

        mGroupAllocPercents.add(allocPerc);
    }

    public final double[] getBucketCounts() { return mBucketCounts; }
    public final double[] getElevatedBucketCounts() { return mElevBucketCounts; }
    public final double[] getNoiseCounts() { return mNoiseCounts; }
    public final double[] getAllocNoiseCounts() { return mAllocNoiseCounts; }
    public double getTotalCount() { return mVarTotal; }
    public double getElevatedCount() { return mElevatedTotal; }
    public double getBackgroundCount() { return mBackgroundTotal; }
    public double getAllocatedCount() { return mAllocTotal; }
    public double getUnallocatedCount() { return mUnallocTotal; }
    public double getAllocNoise() { return mNoiseAllocTotal; }
    public double getMaxNoise() { return mMaxNoiseTotal; }

    public final double[] getUnallocBucketCounts() { return mUnallocBucketCounts; }
    public final double[] getAllocBucketCounts() { return mAllocBucketCounts; }
    public double getAllocPercent() { return mAllocPercent; }
    public double lastAllocPercChange() { return mAllocPercent - mPreviousAllocPerc; }
    public double getUnallocPercent() { return 1 - mAllocPercent; }
    public double getNoisePerc() { return mNoiseAllocTotal/mMaxNoiseTotal; }
    public double getNoiseOfTotal() { return mNoiseAllocTotal/mElevatedTotal; }
    public boolean usingElevatedForAllocation() { return mUseElevatedForAllocation; }

    public void setBucketCounts(final double[] counts)
    {
        mBucketCounts = new double[counts.length];
        copyVector(counts, mBucketCounts);
        mVarTotal = sumVector(counts);
        mMaxNoiseTotal = mVarTotal * MAX_NOISE_ALLOC_PERCENT;
    }

    public void setElevatedBucketCounts(final double[] counts, final double[] noise)
    {
        mElevBucketCounts = new double[counts.length];
        mAllocBucketCounts = new double[counts.length];
        mAllocNoiseCounts = new double[counts.length];
        mUnallocBucketCounts = new double[counts.length];
        mNoiseCounts = new double[counts.length];
        copyVector(counts, mUnallocBucketCounts);
        copyVector(counts, mElevBucketCounts);
        mElevatedTotal = sumVector(counts);
        mBackgroundTotal = mVarTotal - mElevatedTotal;
        mUnallocTotal = mElevatedTotal;

        copyVector(noise, mNoiseCounts);
    }

    public void clearAllocations(boolean useElevatedOnly)
    {
        mBucketGroups.clear();

        // don't force the BG group into the list in case it's never actually allocated to it
        // if(mBackgroundGroup != null)
        //    mBucketGroups.add(mBackgroundGroup);

        mElevBucketGroups.clear();
        mUnallocBuckets.clear();
        mUnallocBuckets.addAll(mElevatedBuckets);

        if(useElevatedOnly)
        {
            copyVector(mElevBucketCounts, mUnallocBucketCounts);
            mUnallocTotal = mElevatedTotal;
        }
        else
        {
            // all allocations from now on will be taken from the full set of counts, not just the elevated ones
            mUseElevatedForAllocation = false;
            copyVector(mBucketCounts, mUnallocBucketCounts);
            mUnallocTotal = mVarTotal;
        }

        initVector(mAllocBucketCounts, 0);
        initVector(mAllocNoiseCounts, 0);
        mGroupAllocPercents.clear();
        mAllocPercent = 0;
        mPreviousAllocPerc = 0;
        mAllocTotal = 0;
        mNoiseAllocTotal = 0;
    }

    public void populateBucketCountSubset(double[] counts, final List<Integer> bucketSubset)
    {
        // extract the counts for the specified subset, leaving the rest zeroed
        initVector(counts, 0);

        for(Integer bucketId : bucketSubset)
        {
            counts[bucketId] = max(mUnallocBucketCounts[bucketId], 0);
        }
    }

    public double getPotentialUnallocCounts(final double[] bucketRatios, final List<Integer> requiredBuckets, final double[] ratioRanges,
            double[] allocCounts)
    {
        return getPotentialCounts(bucketRatios, true, requiredBuckets, ratioRanges, allocCounts);
    }

    public double getPotentialCounts(final double[] bucketRatios, final List<Integer> requiredBuckets,
            final double[] ratioRanges, double[] allocCounts)
    {
        return getPotentialCounts(bucketRatios, false, requiredBuckets, ratioRanges, allocCounts);
    }

    // These next 2 methods are critical logic in how bucket ratios are applied to sample counts, and they work as a pair
    // Firstly a sample's (unallocated) counts are tested against the set of external bucket ratios, making use of both noise in the
    // counts and ratio ranges if applicable. Noise constraints are imposed per bucket, as a proportion of the sample total, and relative
    // to how much of the sample's total this ratios are attempting to allocate

    // Once the potential allocations have been determined, it is logically consistent that they can be applied to the sample without
    // any further constraints or reductions - ie they will fall within remaining unallocated counts and permitted noise
    private double getPotentialCounts(final double[] bucketRatios, boolean useUnallocated, final List<Integer> requiredBuckets,
            final double[] ratioRanges, double[] allocCounts)
    {
        // must cap at the actual sample counts
        // allow to go as high as the elevated probability range
        // if a single bucket required by the ratios has zero unallocated, then all are zeroed
        final double[] sampleCounts = useUnallocated ? mUnallocBucketCounts : (mUseElevatedForAllocation ? mElevBucketCounts : mBucketCounts);
        double countsTotal = mUseElevatedForAllocation ? mElevatedTotal : mVarTotal;
        double unallocTotal = useUnallocated ? mUnallocTotal : countsTotal;

        if(allocCounts != null)
            initVector(allocCounts, 0);

        double allocTotal = 0; // sum of the calculated bucket allocations

        // first extract the remaining unallocated counts per bucket
        double minAlloc = -1;
        double minAllocNoNoise = -1;

        for(Integer bucket : requiredBuckets)
        {
            double unallocCount = max(sampleCounts[bucket], 0);

            if(unallocCount == 0)
            {
                /* this logic is not longer applicable since noise is more carefully managed now

                if(mElevBucketCounts[bucket] > 0)
                {
                    // if any of the required elevated buckets are already fully allocated, no others can be allocated,
                    // regardless of unallocated noise
                    minAlloc = 0;
                    break;
                }
                */

                // otherwise it is valid to use non-negative noise with a zero elevated count
            }

            double noiseCount = useUnallocated ? max(mNoiseCounts[bucket] - mAllocNoiseCounts[bucket], 0) : mNoiseCounts[bucket];

            // allow full allocation of noise, not proportional to allocated counts to make consistent with fitter
            double potentialAlloc = unallocCount + noiseCount;

            // use the low range ratio to give max possible allocation to this bucket
            double ratio = ratioRange(bucketRatios, ratioRanges, bucket, true);

            if(ratio <= 0)
                continue;

            double alloc = potentialAlloc / ratio;
            double allocNoNoise = unallocCount / ratio;

            if(minAlloc == -1 || alloc < minAlloc)
                minAlloc = alloc;

            if(minAllocNoNoise == -1 || allocNoNoise < minAllocNoNoise)
                minAllocNoNoise = allocNoNoise;
        }

        if(minAlloc <= 0)
            return 0;

        double currentNoiseTotal = useUnallocated ? mNoiseAllocTotal : 0;

        // noise is capped at a maximum but also limited by the proportion being allocated to this bucket group
        // eg if maxNoise = 20 but this group is only asking to alloc 50 of 100 total counts, then only 10 of noise can be allocated
        double maxNoiseAllocation = min(mMaxNoiseTotal - currentNoiseTotal, min(minAlloc/countsTotal, 1.0) * mMaxNoiseTotal);

        if(minAllocNoNoise >= mUnallocTotal)
        {
            // special case where no noise is required to allocate remaining counts, so don't over-allocate purely to noise
            minAlloc = minAllocNoNoise;
            maxNoiseAllocation = 0;
        }

        // cap by the actual unallocated before working out allocation
        minAlloc = capValue(minAlloc, 0, unallocTotal + maxNoiseAllocation);

        // if noise needs to be allocated, do this proportionally by ratio amongst the buckets which need it
        double[] bucketWithNoiseRatios = new double[requiredBuckets.size()];

        if(maxNoiseAllocation > 0)
        {
            for (int i = 0; i < requiredBuckets.size(); ++i)
            {
                Integer bucket = requiredBuckets.get(i);

                double ratio = ratioRange(bucketRatios, ratioRanges, bucket, false);
                double allocCount = minAlloc * ratio;

                if (allocCount == 0)
                    continue;

                double noiseCount = useUnallocated ? max(mNoiseCounts[bucket] - mAllocNoiseCounts[bucket], 0) : mNoiseCounts[bucket];

                if (greaterThan(allocCount, sampleCounts[bucket] + noiseCount))
                    allocCount = sampleCounts[bucket] + noiseCount; // shouldn't happen since should have determined above

                if (noiseCount > 0 && sampleCounts[bucket] < allocCount)
                {
                    bucketWithNoiseRatios[i] = ratio;
                }
            }
        }

        double noiseRatioTotal = sumVector(bucketWithNoiseRatios);

        double currentNoiseAlloc = 0;
        for(int i = 0; i < requiredBuckets.size(); ++i)
        {
            Integer bucket = requiredBuckets.get(i);
            double ratio = ratioRange(bucketRatios, ratioRanges, bucket, false);
            double allocCount = minAlloc * ratio;

            if(allocCount == 0)
                continue;

            double noiseCount = useUnallocated ? max(mNoiseCounts[bucket] - mAllocNoiseCounts[bucket], 0) : mNoiseCounts[bucket];

            // allocate remaining noise proportionally by bucket ratio
            double noiseRatio = bucketWithNoiseRatios[i];
            if(noiseRatio > 0)
                noiseCount = min(noiseRatio / noiseRatioTotal * maxNoiseAllocation, noiseCount);

            if(noiseCount + currentNoiseAlloc > maxNoiseAllocation)
            {
                noiseCount = max(maxNoiseAllocation - currentNoiseAlloc, 0);
            }

            if(greaterThan(allocCount, sampleCounts[bucket] + noiseCount))
                allocCount = sampleCounts[bucket] + noiseCount; // check cannot exceed actual + permitted noise

            if(allocCounts != null)
                allocCounts[bucket] = allocCount;

            allocTotal += allocCount;

            if(noiseCount > 0 && sampleCounts[bucket] < allocCount)
            {
                currentNoiseAlloc += allocCount - sampleCounts[bucket];
            }
        }

        // avoid tiny cumulative allocations
        if(allocTotal < 1)
        {
            allocTotal = 0;

            if(allocCounts != null)
                initVector(allocCounts, 0);
        }

        return allocTotal;
    }

    public double allocateBucketCounts(double[] counts, double reqAllocationPercent)
    {
        return allocateBucketCounts(counts, reqAllocationPercent, false);
    }

    public double allocateBucketCounts(double[] counts)
    {
        return allocateBucketCounts(counts, 0, true);
    }

    private double allocateBucketCounts(double[] counts, double reqAllocationPercent, boolean expectFullAllocation)
    {
        double allocatedCount = 0;
        double workingNoiseTotal = mNoiseAllocTotal;
        double refVarTotal = mUseElevatedForAllocation ? mElevatedTotal : mVarTotal;

        if(!expectFullAllocation)
        {
            double reductionFactor = 1;

            // do a preliminary check that this allocation will actually achieve the required percentage count
            for (int i = 0; i < counts.length; ++i)
            {
                if (counts[i] == 0)
                    continue;

                double unallocNoise = max(mNoiseCounts[i] - mAllocNoiseCounts[i], 0);

                // the total allocated to noise across all buckets cannot exceed a factor of the sample total
                if (unallocNoise + workingNoiseTotal > mMaxNoiseTotal)
                {
                    unallocNoise = max(mMaxNoiseTotal - workingNoiseTotal, 0);
                }

                double allocCount = min(mUnallocBucketCounts[i] + unallocNoise, counts[i]);

                reductionFactor = min(reductionFactor, allocCount / counts[i]);

                // if even a single bucket restricts a non-zero count being allocated, fail the whole allocation
                if (allocCount == 0)
                    return 0;

                allocatedCount += allocCount;

                if (unallocNoise > 0 && mUnallocBucketCounts[i] < allocCount)
                {
                    workingNoiseTotal += allocCount - mUnallocBucketCounts[i];
                }
            }

            if (lessThan(reductionFactor, 1))
            {
                // reduce counts proportionally since they are generated from ratios
                allocatedCount *= reductionFactor;
                vectorMultiply(counts, reductionFactor);
            }

            // now only factors in the allocation to remaining actual counts
            if (reqAllocationPercent > 0 && allocatedCount < reqAllocationPercent * refVarTotal)
                return allocatedCount;
        }

        // allow allocation for the caller to go up to the unallocated counts plus the permitted noise range, which
        // will only be allocated once (when count > unallocated)
        // modify the caller's count values if limited in any way
        // internal allocation counts and total are limited to actuals
        workingNoiseTotal = mNoiseAllocTotal;
        allocatedCount = 0;

        for(int i = 0; i < counts.length; ++i)
        {
            if(counts[i] == 0)
                continue;

            if (mUnallocBucketCounts[i] < counts[i])
            {
                // a greater allocation has been requested than is available, therefore look to cover with unallocated noise

                // first limit by the bucket
                double maxUnallocBucketNoise = max(mNoiseCounts[i] - mAllocNoiseCounts[i], 0);

                // then by the max per sample
                if(maxUnallocBucketNoise + workingNoiseTotal > mMaxNoiseTotal)
                {
                    maxUnallocBucketNoise = max(mMaxNoiseTotal - workingNoiseTotal, 0);
                }

                if(mUnallocBucketCounts[i] == 0 && maxUnallocBucketNoise == 0)
                {
                    // neither unallocated actuals counts nor unallocated noise can be used
                    if(expectFullAllocation & counts[i] > 0.01)
                    {
                        LOGGER.warn("sample({}) cannot alloc bucket({}) - no unalloc count or noise", Id, i);
                    }

                    counts[i] = 0;
                    continue;
                }

                // fully allocated now
                double actualAlloc = mUnallocBucketCounts[i];

                if(actualAlloc > 0)
                {
                    mAllocBucketCounts[i] += actualAlloc;
                    mUnallocBucketCounts[i] = 0;
                    removeUnallocBucket(i);
                }

                double requiredNoiseAlloc = counts[i] - actualAlloc;
                double actualNoiseAlloc = min(requiredNoiseAlloc, maxUnallocBucketNoise);

                if(Doubles.lessThan(actualAlloc + actualNoiseAlloc, counts[i]))
                {
                    // may be an issue if more was requested than could be allocated
                    if(expectFullAllocation && greaterThan(actualAlloc + actualNoiseAlloc, counts[i]))
                    {
                        LOGGER.warn(String.format("sample(%s) reducing bucket(%d) since requested(%.3f) < unalloc(%.3f) + noise(%.3f)",
                                Id, i, counts[i], actualAlloc, actualNoiseAlloc));
                    }

                    counts[i] = actualAlloc + actualNoiseAlloc;
                }

                mNoiseAllocTotal += actualNoiseAlloc;
                mAllocNoiseCounts[i] += actualNoiseAlloc;
                workingNoiseTotal += actualNoiseAlloc;
            }
            else
            {
                mUnallocBucketCounts[i] -= counts[i];
                mAllocBucketCounts[i] += counts[i];
            }

            allocatedCount += counts[i];
        }

        mAllocTotal = capValue(sumVector(mAllocBucketCounts), 0, refVarTotal);
        mUnallocTotal = refVarTotal - mAllocTotal;

        mPreviousAllocPerc = mAllocPercent;
        mAllocPercent = capValue(mAllocTotal/refVarTotal, 0, 1);

        return allocatedCount;
    }

    public void restoreCounts(final double[] allocCounts, final double[] noiseCounts)
    {
        double refVarTotal = mUseElevatedForAllocation ? mElevatedTotal : mVarTotal;
        copyVector(allocCounts, mAllocBucketCounts);
        copyVector(noiseCounts, mAllocNoiseCounts);
        mAllocTotal = 0;
        mNoiseAllocTotal = 0;

        for(int i = 0; i < allocCounts.length; ++i)
        {
            double bucketCount = mUseElevatedForAllocation ? mElevBucketCounts[i] : mBucketCounts[i];
            mUnallocBucketCounts[i] = bucketCount - mAllocBucketCounts[i];
            mAllocTotal += mAllocBucketCounts[i];
            mNoiseAllocTotal += mAllocNoiseCounts[i];

            checkAddUnallocBucket(i);
        }

        mUnallocTotal = refVarTotal - mAllocTotal;
        mAllocPercent = capValue(mAllocTotal/refVarTotal, 0, 1);
    }

    public double reduceAllocCounts(final double[] counts)
    {
        double refVarTotal = mUseElevatedForAllocation ? mElevatedTotal : mVarTotal;
        double reductionTotal = 0;

        for(int i = 0; i < counts.length; ++i)
        {
            if(counts[i] >= 0)
                continue;

            // take from allocated noise first of all if allocated counts are at their limit
            double countReduction = -counts[i];

            if(mUnallocBucketCounts[i] == 0)
            {
                double noiseReduction = 0;
                double noiseAlloc = mAllocNoiseCounts[i];

                if(noiseAlloc >= countReduction)
                {
                    noiseReduction = countReduction;
                    countReduction = 0;
                }
                else
                {
                    noiseReduction = noiseAlloc;
                    countReduction = countReduction - noiseAlloc;
                }

                if(noiseReduction > 0)
                {
                    mAllocNoiseCounts[i] -= noiseReduction;
                    mNoiseAllocTotal -= noiseReduction;
                    reductionTotal += noiseReduction;
                }
            }

            if(countReduction > 0)
            {
                mUnallocBucketCounts[i] += countReduction;
                mAllocBucketCounts[i] -= countReduction;
                reductionTotal += countReduction;
                checkAddUnallocBucket(i);
            }
        }

        mAllocTotal = max(mAllocTotal - reductionTotal, 0);
        mUnallocTotal = refVarTotal - mAllocTotal;

        mPreviousAllocPerc = mAllocPercent;
        mAllocPercent = capValue(mAllocTotal/refVarTotal, 0, 1);

        return -reductionTotal;
    }


    private void removeUnallocBucket(Integer bucketId)
    {
        mUnallocBuckets.remove(bucketId);
    }

    private void checkAddUnallocBucket(int bucketId)
    {
        if(mElevatedBuckets.contains(bucketId) && !mUnallocBuckets.contains(bucketId))
        {
            if(mUnallocBucketCounts[bucketId] > 0)
                mUnallocBuckets.add(bucketId);
        }
    }

}
