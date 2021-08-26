package com.hartwig.hmftools.sigs.fitter;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.sigs.DataUtils.doubleToStr;
import static com.hartwig.hmftools.common.sigs.SigUtils.calcResiduals;
import static com.hartwig.hmftools.common.utils.VectorUtils.copyVector;
import static com.hartwig.hmftools.common.utils.VectorUtils.initVector;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;
import static com.hartwig.hmftools.sigs.buckets.BucketGroup.ratioRange;
import static com.hartwig.hmftools.sigs.common.CommonUtils.SIG_LOGGER;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sigs.SigResiduals;
import com.hartwig.hmftools.common.utils.Doubles;

public class ConstrainedFitter
{
    private final double[] mCounts; // sample counts, optionally including noise
    private final double[] mAllocCounts; // counts allocated from signature contributions

    // private int mSampleId;

    private double[] mContribs;
    private double mContribTotal;

    private final List<double[]> mRatiosCollection; // the signature bucket ratios
    private final List<String> mSigNames;
    private final List<Integer> mZeroedSigs;

    private final int mBucketCount;
    private int mSigCount;
    private double mCountsTotal; // total of the actual variant counts

    // parameters
    private double mTargetAllocPercent; // target total allocation to exit the fit
    private double mMinAllocation;
    private final List<Integer> mRequiredSigs;
    private SigResiduals mResiduals;

    private boolean mIsValid;
    private boolean mLogVerbose;

    public ConstrainedFitter(int bucketCount)
    {
        mBucketCount = bucketCount;

        mCounts = new double[mBucketCount];
        mAllocCounts = new double[mBucketCount];

        mLogVerbose = false;
        mTargetAllocPercent = 0.001;

        mRatiosCollection = Lists.newArrayList();
        mSigNames = Lists.newArrayList();
        mZeroedSigs = Lists.newArrayList();
        mRequiredSigs = Lists.newArrayList();
        mSigCount = 0;
        mMinAllocation = 0;
        mResiduals = new SigResiduals();

        clearData();

        mIsValid = true;
    }

    private void clearData()
    {
        mContribs = null;
        mRatiosCollection.clear();
        mSigNames.clear();
        mSigCount = 0;
        mZeroedSigs.clear();
        mRequiredSigs.clear();

        mCountsTotal = 0;
        mContribTotal = 0;
        mResiduals = new SigResiduals();

        initVector(mAllocCounts, 0);
    }

    public void setParameters(double minAllocation, final List<Integer> requiredSigs, boolean logVerbose, double targetAllocPercent)
    {
        mMinAllocation = minAllocation;

        if(requiredSigs != null)
            mRequiredSigs.addAll(requiredSigs);

        mLogVerbose = logVerbose;
        mTargetAllocPercent = targetAllocPercent;
    }

    public final double[] getContribs() { return mContribs; }
    public double getAllocPerc() { return mContribTotal / mCountsTotal; }
    public boolean isValid() { return mIsValid; }

    public boolean fitCounts(final List<double[]> ratiosCollection, final double[] counts, final List<String> sigNames)
    {
        mIsValid = true;
        clearData();

        mCountsTotal = sumVector(counts);

        if(ratiosCollection.isEmpty() || mCountsTotal <= 0 || counts.length != mBucketCount)
        {
            mIsValid = false;
            return false;
        }
        else if(ratiosCollection.stream().anyMatch(x -> x.length != mBucketCount))
        {
            mIsValid = false;
            return false;
        }

        copyVector(counts, mCounts);

        mRatiosCollection.addAll(ratiosCollection);
        mSigCount = mRatiosCollection.size();

        if(sigNames != null && sigNames.size() == mSigCount)
            mSigNames.addAll(sigNames);

        mContribs = new double[mSigCount];

        initialAllocation();
        logStats();

        // reduce each exhausted bucket by the most cost-effective manner
        reduceExhaustedBuckets();
        logStats();

        // finally attempt re-assignment for each signature up to the permitted levels
        finalAllocation();
        logStats();

        return true;
    }

    private void allocateByCost()
    {
        // allocate to the smallest ratios first up to their maximum
    }

    private void initialAllocation()
    {
        // add each signature up to its maximum

        for(int s = 0; s < mSigCount; ++s)
        {
            final double[] sigRatios = mRatiosCollection.get(s);

            double minAlloc = -1;

            for(int b = 0; b < mBucketCount; ++b)
            {
                double count = mCounts[b];

                if(sigRatios[b] == 0)
                    continue;

                if(count == 0)
                {
                    minAlloc = 0;
                    break;
                }

                double maxAlloc = count / sigRatios[b];

                minAlloc = minAlloc < 0 ? maxAlloc : min(minAlloc, maxAlloc);
            }

            if(minAlloc == 0)
                continue;

            if(mLogVerbose)
            {
                SIG_LOGGER.debug("sig({}) initial allocation({})", getSigName(s), doubleToStr(minAlloc));
            }

            applySignatureContribution(s, minAlloc);
        }
    }

    private static final int BUCKET_ID = 0;
    private static final int BUCKET_EXCESS = 1;

    private void reduceExhaustedBuckets()
    {
        final List<int[]> exhaustedBucketData = Lists.newArrayList();

        for (int b = 0; b < mBucketCount; ++b)
        {
            if (mAllocCounts[b] <= mCounts[b])
                continue;

            int excess = (int)(mAllocCounts[b] - mCounts[b]);

            int index = 0;
            while(index < exhaustedBucketData.size())
            {
                if(excess > exhaustedBucketData.get(index)[BUCKET_EXCESS])
                    break;


                ++index;
            }

            exhaustedBucketData.add(index, new int[] {b, excess});
        }

        if(exhaustedBucketData.isEmpty())
            return;

        SIG_LOGGER.debug("found {} excess buckets", exhaustedBucketData.size());

        for(final int[] bucketData : exhaustedBucketData)
        {
            int bucket = bucketData[BUCKET_ID];

            double requiredReduction = mAllocCounts[bucket] - mCounts[bucket];

            if(requiredReduction <= 0)
                continue;

            // from all contributing signatures, find the relative reduction amounts for each
            double[] requiredReductions = new double[mSigCount];
            double[] sigReductionCosts = new double[mSigCount];
            double totalCost = 0;

            for(int s = 0; s < mSigCount; ++s)
            {
                final double[] sigRatios = mRatiosCollection.get(s);

                if(sigRatios[bucket] == 0)
                    continue;

                double cost = 1 / sigRatios[bucket];
                double ratioCost = sigRatios[bucket];
                sigReductionCosts[s] = ratioCost;
                totalCost += ratioCost;
            }

            for(int s = 0; s < mSigCount; ++s)
            {
                if(sigReductionCosts[s] == 0)
                    continue;

                final double[] sigRatios = mRatiosCollection.get(s);
                double sigBucketReduction = requiredReduction * sigReductionCosts[s] / totalCost;
                double sigReduction = sigBucketReduction / sigRatios[bucket];

                applySignatureContribution(s, -sigReduction);
            }

            if(!Doubles.equal(mAllocCounts[bucket], mCounts[bucket]))
            {
                SIG_LOGGER.error("bucket({}) reduction failed", bucket);
            }
        }
    }

    private void finalAllocation()
    {
        // add each signature up to its maximum
        for(int s = 0; s < mSigCount; ++s)
        {
            final double[] sigRatios = mRatiosCollection.get(s);

            double minAlloc = -1;

            for(int b = 0; b < mBucketCount; ++b)
            {
                double unallocCount = mCounts[b] - mAllocCounts[b];

                if(sigRatios[b] == 0)
                    continue;

                if(unallocCount <= 0)
                {
                    minAlloc = 0;
                    break;
                }

                double maxAlloc = unallocCount / sigRatios[b];

                minAlloc = minAlloc < 0 ? maxAlloc : min(minAlloc, maxAlloc);
            }

            if(minAlloc == 0)
                continue;

            if(mLogVerbose)
            {
                SIG_LOGGER.debug("sig({}) final allocation({})", getSigName(s), doubleToStr(minAlloc));
            }

            applySignatureContribution(s, minAlloc);
        }
    }

    private void applySignatureContribution(int sig, double allocation)
    {
        final double[] sigRatios = mRatiosCollection.get(sig);

        for(int b = 0; b < mBucketCount; ++b)
        {
            mAllocCounts[b] += sigRatios[b] * allocation;

            if(mAllocCounts[b] < 0)
            {
                SIG_LOGGER.error("sig({}) bucket({}) count reduced below zero for allocation({})",
                        getSigName(sig), b, doubleToStr(allocation));
                mIsValid = false;
                return;
            }
        }

        mContribTotal += allocation;
        mContribs[sig] += allocation;
    }

    private void logStats()
    {
        mResiduals = calcResiduals(mCounts, mAllocCounts, mCountsTotal);

        if(!mLogVerbose)
            return;

        SIG_LOGGER.debug(String.format("totalCount(%s) allocated(%s perc=%.3f) res(%s perc=%.3f) underAlloc(%s perc=%.3f) excess(%s perc=%.3f)",
                doubleToStr(mCountsTotal), doubleToStr(mContribTotal), mContribTotal / mCountsTotal,
                doubleToStr(mResiduals.Total), mResiduals.Percent, doubleToStr(mResiduals.unallocated()), mResiduals.unallocated() / mCountsTotal,
                doubleToStr(mResiduals.Excess), mResiduals.Excess / mCountsTotal));
    }

    private String getSigName(int sig) { return !mSigNames.isEmpty() ? mSigNames.get(sig) : String.valueOf(sig); }

}
