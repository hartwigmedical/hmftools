package com.hartwig.hmftools.sig_analyser.buckets;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.sig_analyser.buckets.BaConfig.MAX_NOISE_ALLOC_PERCENT;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.capValue;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.copyVector;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.greaterThan;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.initVector;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.lessThan;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.sumVector;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.vectorMultiply;
import static com.hartwig.hmftools.sig_analyser.buckets.BucketGroup.ratioRange;

import java.util.List;

import com.google.common.collect.Lists;

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
    private double[] mAllocNoiseCounts; // counts allocated to noise for each bucket
    private double[] mAllocBucketCounts; // allocated counts for elevated buckets (excluding noise)
    private double[] mUnallocBucketCounts; // unallocated counts for elevated buckets
    private double mAllocPercent;
    private double mPreviousAllocPerc; // to allow tracking of changes
    private double mVarTotal;
    private double mBackgroundTotal;
    private double mElevatedTotal;
    private double mAllocTotal;
    private double mUnallocTotal;
    private double mNoiseAllocTotal;
    private boolean mAllocOfElevated; // any calculated allocations will by from elevated, not total count
    private double mMaxNoiseTotal;

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
        mAllocOfElevated = true;
        mPreviousAllocPerc = 0;
        mAllocTotal = 0;
        mElevatedTotal = 0;
        mBackgroundTotal = 0;
        mUnallocTotal = 0;
        mNoiseAllocTotal = 0;
        mVarTotal = 0;
        mMaxNoiseTotal = 0;
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
    public final double[] getCountRanges() { return mNoiseCounts; }
    public final double[] getAllocNoiseCounts() { return mAllocNoiseCounts; }
    public double getTotalCount() { return mVarTotal; }
    public double getElevatedCount() { return mElevatedTotal; }
    public double getBackgroundCount() { return mBackgroundTotal; }
    public double getAllocatedCount() { return mAllocTotal; }
    public double getUnallocatedCount() { return mUnallocTotal; }
    public double getAllocNoise() { return mNoiseAllocTotal; }
    public double getMaxNoise() { return mMaxNoiseTotal; }

    public final double[] getUnallocBucketCounts() { return mUnallocBucketCounts; }
    public double getAllocPercent() { return mAllocPercent; }
    public double lastAllocPercChange() { return mAllocPercent - mPreviousAllocPerc; }
    public double getUnallocPercent() { return 1 - mAllocPercent; }
    public double getNoisePerc() { return mNoiseAllocTotal/mMaxNoiseTotal; }
    public double getNoiseOfTotal() { return mNoiseAllocTotal/mElevatedTotal; }

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

        if(mBackgroundGroup != null)
            mBucketGroups.add(mBackgroundGroup);

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
            mAllocOfElevated = false;
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

    public double[] getPotentialUnallocCounts(final double[] bucketRatios, final List<Integer> requiredBuckets, final double[] ratioRanges)
    {
        return getPotentialAllocation(bucketRatios, true, requiredBuckets, ratioRanges);
    }

    public double[] getPotentialElevCounts(final double[] bucketRatios, final List<Integer> requiredBuckets, final double[] ratioRanges)
    {
        return getPotentialAllocation(bucketRatios, false, requiredBuckets, ratioRanges);
    }

    private double[] getPotentialAllocation(final double[] bucketRatios, boolean useUnallocated, final List<Integer> requiredBuckets,
            final double[] ratioRanges)
    {
        // must cap at the actual sample counts
        // allow to go as high as the elevated probability range
        // if a single bucket required by the ratios has zero unallocated, then all are zeroed
        final double[] sampleCounts = useUnallocated ? mUnallocBucketCounts : mElevBucketCounts;
        double countsTotal = useUnallocated ? mUnallocTotal : (mAllocOfElevated ? mElevatedTotal : mVarTotal);

        // first extract the remaining unallocated counts per bucket
        double minAlloc = 0;
        double minAllocNoNoise = 0;

        for(Integer bucket : requiredBuckets)
        {
            double unallocCount = max(sampleCounts[bucket], 0);

            if(unallocCount == 0)
            {
                if(mElevBucketCounts[bucket] > 0)
                {
                    // if any of the required elevated buckets are already fully allocated, no others can be allocated,
                    // regardless of unallocated noise
                    minAlloc = 0;
                    break;
                }

                // otherwise it is valid to use non-negative noise with a zero elevated count
            }

            double noiseCount = useUnallocated ? max(mNoiseCounts[bucket] - mAllocNoiseCounts[bucket], 0) : mNoiseCounts[bucket];

            // allow full allocation of noise, not proportional to allocated counts to make consistent with fitter
            double potentialAlloc = unallocCount + noiseCount;

            // use the low range ratio to give max possible allocation to this bucket
            double ratio = bucketRatios[bucket] + ratioRange(ratioRanges, bucket, true);
            double alloc = potentialAlloc / ratio;
            double allocNoNoise = unallocCount / ratio;

            if(minAlloc == 0 || alloc < minAlloc)
                minAlloc = alloc;

            if(minAllocNoNoise == 0 || allocNoNoise < minAllocNoNoise)
                minAllocNoNoise = allocNoNoise;
        }

        double[] bestAllocCounts = new double[bucketRatios.length];

        if(minAlloc == 0)
            return bestAllocCounts;

        double currentNoiseTotal = useUnallocated ? mNoiseAllocTotal : 0;

        // noise is capped at a maximum but also limited by the proportion being allocated to this bucket group
        // eg if maxNoise = 20 but this group is only asking to alloc 50 of 100 total counts, then only 10 of noise can be allocated
        double maxNoiseAllocation = min(mMaxNoiseTotal - currentNoiseTotal, minAllocNoNoise/countsTotal * mMaxNoiseTotal);

        // cap by the actual unallocated before working out allocation
        minAlloc = capValue(minAlloc, 0, countsTotal + maxNoiseAllocation);

        // if noise needs to be allocated, do this proportionally by ratio amongst the buckets which need it
        List<Double> bucketWithNoiseRatios = Lists.newArrayList();

        for(int i = 0; i < requiredBuckets.size(); ++i)
        {
            Integer bucket = requiredBuckets.get(i);
            bucketWithNoiseRatios.add(0.0);

            double ratio = bucketRatios[bucket] + ratioRange(ratioRanges, bucket, false);
            double allocCount = minAlloc * ratio;

            if(allocCount == 0)
                continue;

            double noiseCount = useUnallocated ? max(mNoiseCounts[bucket] - mAllocNoiseCounts[bucket], 0) : mNoiseCounts[bucket];

            if(greaterThan(allocCount, sampleCounts[bucket] + noiseCount))
                allocCount = sampleCounts[bucket] + noiseCount; // shouldn't happen since should have determined above

            if(noiseCount > 0 && sampleCounts[bucket] < allocCount)
            {
                bucketWithNoiseRatios.set(i, ratio);
            }
        }

        double noiseRatioTotal = bucketWithNoiseRatios.stream().mapToDouble(x -> x).sum();

        double noiseAlloc = 0;
        for(int i = 0; i < requiredBuckets.size(); ++i)
        {
            Integer bucket = requiredBuckets.get(i);
            double ratio = bucketRatios[bucket] + ratioRange(ratioRanges, bucket, false);
            double allocCount = minAlloc * ratio;

            if(allocCount == 0)
                continue;

            double noiseCount = useUnallocated ? max(mNoiseCounts[bucket] - mAllocNoiseCounts[bucket], 0) : mNoiseCounts[bucket];

            // allocate remaining noise proportionally by bucket ratio
            double noiseRatio = bucketWithNoiseRatios.get(i);
            if(noiseRatio > 0)
                noiseCount = noiseRatio / noiseRatioTotal * maxNoiseAllocation;

            if(noiseCount + noiseAlloc > maxNoiseAllocation)
            {
                noiseCount = max(maxNoiseAllocation - noiseAlloc, 0);
            }

            if(greaterThan(allocCount, sampleCounts[bucket] + noiseCount))
                allocCount = sampleCounts[bucket] + noiseCount; // shouldn't happen since should have determined above

            bestAllocCounts[bucket] = allocCount;

            if(noiseCount > 0 && sampleCounts[bucket] < allocCount)
            {
                noiseAlloc += allocCount - sampleCounts[bucket];
            }
        }

        return bestAllocCounts;
    }

    public double allocateBucketCounts(double[] counts, double reqAllocationPercent)
    {
        double allocatedCount = 0;
        double noiseTotal = mNoiseAllocTotal;
        double refVarTotal = mAllocOfElevated ? mElevatedTotal : mVarTotal;
        double reductionFactor = 1;

        // do a preliminary check that this allocation will actually achieve the required percentage count
        for(int i = 0; i < counts.length; ++i)
        {
            if(counts[i] == 0)
                continue;

            double unallocNoise = max(mNoiseCounts[i] - mAllocNoiseCounts[i], 0);

            // the total allocated to noise across all buckets cannot exceed a factor of the sample total
            if(unallocNoise + noiseTotal > mVarTotal * MAX_NOISE_ALLOC_PERCENT)
            {
                unallocNoise = max(mMaxNoiseTotal - noiseTotal, 0);
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
        }

        if(lessThan(reductionFactor, 1))
        {
            // assumption is that the counts are in proportion (eg if generated off
            allocatedCount *= reductionFactor;
            vectorMultiply(counts, reductionFactor);
        }

        // now only factors in the allocation to remaining actual counts
        if(reqAllocationPercent > 0 && allocatedCount < reqAllocationPercent * refVarTotal)
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
                double unallocNoise = max(mNoiseCounts[i] - mAllocNoiseCounts[i], 0);

                if(unallocNoise + noiseTotal > mMaxNoiseTotal)
                {
                    unallocNoise = max(mMaxNoiseTotal - noiseTotal, 0);
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
        }

        mAllocTotal = capValue(sumVector(mAllocBucketCounts), 0, refVarTotal);
        mUnallocTotal = refVarTotal - mAllocTotal;

        mPreviousAllocPerc = mAllocPercent;
        mAllocPercent = capValue(mAllocTotal/refVarTotal, 0, 1);

        return allocatedCount;
    }

    private void removeUnallocBucket(Integer bucketId)
    {
        mUnallocBuckets.remove(bucketId);
    }

}
