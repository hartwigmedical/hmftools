package com.hartwig.hmftools.data_analyser.calcs;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.numeric.Doubles.lessThan;
import static com.hartwig.hmftools.common.numeric.Doubles.greaterThan;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.capValue;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.copyVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.doublesEqual;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.initVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sizeToStr;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sumVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sumVectors;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.data_analyser.types.NmfMatrix;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SigContributionOptimiser
{
    private int mSampleId;
    private boolean mIsValid;

    private double[] mRawCounts; // actual sample counts
    private double[] mCountsNoise; // noise around the counts
    private List<double[]> mRatiosCollection; // the signature bucket ratios

    private double[] mInitContribs;
    private double mMinContribPercent; // each sig's final contrib must be above this level
    private double mTargetAllocPercent; // target total allocation to exit the fit
    private int mTargetSig; // to check if a specific sig remains above the require contribution percent
    private int mRequiredSig; // keep this specific sig even if it falls below the require contribution percent
    private List<Integer> mZeroedSigs; // the signature bucket ratios

    private double[] mCounts; // sample counts, optionally including noise
    private double[] mCurrentCounts;
    private double[] mContribs;
    private double mContribTotal;

    private NmfMatrix mSigs;
    private NmfMatrix mSigCostBasis;

    private int mBucketCount;
    private int mSigCount;
    private double mCountsTotal; // total of the counts including any applied noise
    private double mVarCount; // total of the actual variant counts
    private double mResiduals;
    private double mCurrentAllocPerc;
    private double mInitAllocPerc;
    private double mMinContribChange;
    private boolean mHasLowContribSigs;
    private boolean mIsFullyAllocated;

    private boolean mLogVerbose;
    private boolean mApplyNoise;

    // constants and control config
    private static double MIN_COUNT_CHG_PERC = 0.001;
    private static double MIN_COUNT_CHG = 1;

    private static final Logger LOGGER = LogManager.getLogger(SigContributionOptimiser.class);

    public SigContributionOptimiser(int bucketCount, boolean logVerbose, double minRequiredPerc, double targetAllocPercent, boolean applyNoise)
    {
        mBucketCount = bucketCount;

        mLogVerbose = logVerbose;
        mMinContribPercent = minRequiredPerc;
        mTargetAllocPercent = targetAllocPercent;
        mApplyNoise = applyNoise;

        mTargetSig = -1;
        mRequiredSig = -1;
        mZeroedSigs = Lists.newArrayList();
        mRatiosCollection = Lists.newArrayList();

        mCurrentCounts = new double[mBucketCount];
        mCounts = new double[mBucketCount];
    }

    public void initialise(int sampleId, final double[] counts, final double[] countsNoise, final List<double[]> ratiosCollection, final double[] contribs)
    {
        mSampleId = sampleId;

        mSigCount = ratiosCollection.size();

        mRawCounts = counts;
        mCountsNoise = countsNoise;

        copyVector(counts, mCounts);

        mTargetSig = -1;
        mRequiredSig = -1;
        mZeroedSigs.clear();

        sumVectors(mCountsNoise, mCounts);

        mRatiosCollection.clear();
        mRatiosCollection.addAll(ratiosCollection);

        mContribs = new double[mSigCount];
        mInitContribs = new double[mSigCount];
        copyVector(contribs, mInitContribs);
        mContribTotal = 0;

        mVarCount = sumVector(mRawCounts);
        mCountsTotal = sumVector(mCounts);
        mResiduals = 0;
        mCurrentAllocPerc = 0;
        mInitAllocPerc = 0;
        mHasLowContribSigs = false;
        mIsFullyAllocated = false;
        mMinContribPercent = 0;
        mTargetAllocPercent = 0;

        mMinContribChange = max(mVarCount * MIN_COUNT_CHG_PERC, MIN_COUNT_CHG);

        for(int b = 0; b < mBucketCount; ++b)
        {
            mCurrentCounts[b] = 0;
        }

        // extract sigs their cost basis (just the inverse)
        mSigs = new NmfMatrix(mBucketCount, mSigCount);
        mSigCostBasis = new NmfMatrix(mBucketCount, mSigCount);

        if(mContribs.length != ratiosCollection.size())
        {
            mIsValid = false;
            return;
        }

        double[][] sigData = mSigs.getData();
        double[][] cbData = mSigCostBasis.getData();

        for(int sig = 0; sig < mSigCount; ++sig)
        {
            double[] sigRatios = ratiosCollection.get(sig);

            if(sigRatios.length != mBucketCount)
            {
                mIsValid = false;
                return;
            }

            for(int b = 0; b < mBucketCount; ++b)
            {
                sigData[b][sig] = sigRatios[b];

                if(sigRatios[b] > 0)
                    cbData[b][sig] = 1 / sigRatios[b];
            }

            mContribs[sig] = 0;
        }

        mIsValid = true;
    }

    public final double[] getFittedCounts() { return mCurrentCounts; }
    public final double[] getContribs() { return mContribs; }
    public final double getContribTotal() { return mContribTotal; }
    public double getResiduals() { return mResiduals; }
    public double getAllocPerc() { return mCurrentAllocPerc; }
    public boolean isValid() { return mIsValid; }
    public boolean isFullyAllocated() { return mIsFullyAllocated; }

    public void setTargetSig(int sig) { mTargetSig = sig; }
    public void setRequiredSig(int sig) { mRequiredSig = sig; }
    public void setLogVerbose(boolean toggle) { mLogVerbose = toggle; }
    public void setApplyNoise(boolean toggle)
    {
        mApplyNoise = toggle;

        if(!mApplyNoise)
        {
            copyVector(mRawCounts, mCounts);
        }

        mCountsTotal = sumVector(mCounts);
    }

    private void applyInitialContributions()
    {
        if(mRequiredSig >= 0)
            calcSigContribution(mRequiredSig);

        calcAllContributions();

        copyVector(mContribs, mInitContribs);
        double initContrib = sumVector(mInitContribs);
        mInitAllocPerc = initContrib/mVarCount;

        for(int s = 0; s < mSigCount; ++s)
        {
            if(s == mRequiredSig)
                continue;

            if(mInitContribs[s]/mVarCount < mMinContribPercent)
                mHasLowContribSigs = true;
        }
    }

    public boolean fitToSample(double targetAllocPercent, double minSigPercent)
    {
        if (!mIsValid)
            return false;

        mMinContribPercent = minSigPercent;
        mTargetAllocPercent = targetAllocPercent;

        // first apply all the sigs' allocations in turn to exhaust as much of the mCounts as possible
        applyInitialContributions();
        logStats();

        if (!mIsValid)
            return false;

        if(mIsFullyAllocated)
        {
            clearLowContribSigs();
            return true;
        }

        boolean foundImprovements = fitToSample();

        while(foundImprovements || mHasLowContribSigs)
        {
            // strip out the worst contributor if below teh required threshold and try again
            List<Integer> sortedContribIndices = getSortedVectorIndices(mContribs, true);

            for (Integer s : sortedContribIndices)
            {
                if(mContribs[s] == 0)
                    continue;

                if(s == mRequiredSig)
                    continue;

                if(mContribs[s]/mVarCount < mMinContribPercent)
                {
                    zeroSigContrib(s);

                    if(s == mTargetSig)
                    {
                        // LOGGER.debug("sample({}) target sig({}) too low", mSampleId, s);
                        return true;
                    }

                    //if(mLogVerbose)
                    //    LOGGER.debug("sample({}) removed low-percent sig({})", mSampleId, s);

                    calcAllContributions();
                    break; // only remove the worst
                }
            }

            // otherwise try again
            logStats();
            foundImprovements = fitToSample();

            if(mIsFullyAllocated)
            {
                clearLowContribSigs();
                return true;
            }
        }

        return mIsValid;
    }

    private boolean fitToSample()
    {
        // continue optimising until there has been a run through all exhausted buckets with no improvements
        int iterations = 0;
        int maxIterations = 50;

        boolean recalcRequired = false;

        while(iterations < maxIterations)
        {
            List<Integer> exhaustedBuckets = getExhaustedBuckets();

            if(exhaustedBuckets.isEmpty())
            {
                calcAllContributions();
                exhaustedBuckets = getExhaustedBuckets();

                if(exhaustedBuckets.isEmpty())
                    return false;
            }

            List<Integer> sortedContribIndices = getSortedVectorIndices(mContribs, false);

            boolean foundTransfers = false;

            for (Integer eb : exhaustedBuckets)
            {
                for (Integer sig : sortedContribIndices)
                {
                    foundTransfers = reduceSigContribution(sig, eb);

                    if(!mIsValid)
                        return false;

                    if (foundTransfers)
                    {
                        // if an contributions have changed, need to break and begin the process again
                        recalcRequired = true;
                        logStats();
                        break;
                    }
                }

                if (foundTransfers)
                    break;
            }

            if(mIsFullyAllocated)
            {
                if(mLogVerbose)
                    LOGGER.debug(String.format("sample(%d) allocPerc(%.3f) exceeds target", mSampleId, mCurrentAllocPerc));

                break;
            }

            if(!foundTransfers) // no point continuing to look
                break;

            ++iterations;
        }

        return recalcRequired;
    }

    private boolean reduceSigContribution(int rs, int eb)
    {
        /* find the optimal way to reduce each a signature in an exhausted bucket and gain in another signature:
            1. Find all contributing sigs, their contrib and their cost basis
            2. In every other bucket, with this contributiong temporarily removed, find the maximum gain
            3. Apply this maximum gain across all buckets for that sig
            4. Use these to calculate a gain ratio
            5. Select and apply the highest gain ratio
         */

        double[][] sigData = mSigs.getData();
        double[][] cbData = mSigCostBasis.getData();

        if(cbData[eb][rs] == 0)
            return false;

        double[] lessSig1Counts = new double[mBucketCount];

        // calc the minimum the contribution loss in the exhausted bucket for the sig being reduced
        double sig1EBContribLoss = mMinContribChange * cbData[eb][rs]; // mMinContribChange

        sig1EBContribLoss = min(sig1EBContribLoss, mContribs[rs]); // cannot reduce past zero

//        if (sig1EBContribLoss < mMinContribChange)
//            return false;

        double cumulativeSig1ContribLoss = 0;
        double totalSig1ContribLoss = 0;

        while (mContribs[rs] > 0)
        {
            // cache the best option for transfer
            int topSig2 = 0;
            double topSig2ContribGain = 0;

            copyVector(mCurrentCounts, lessSig1Counts);

            cumulativeSig1ContribLoss += sig1EBContribLoss;

            if(mContribs[rs] - cumulativeSig1ContribLoss < 0) // cannot reduce further
                break;

            // remove this sig across the board from a 1-lot contrib to this exhausted bucket
            for (int b = 0; b < mBucketCount; ++b)
            {
                lessSig1Counts[b] -= cumulativeSig1ContribLoss * sigData[b][rs];
            }

            // now look at the potential gain to all the sigs in each bucket, having had sig1 removed
            for (int s2 = 0; s2 < mSigCount; ++s2)
            {
                if (rs == s2)
                    continue;

                // want to find the min allocation for this sig2 from all buckets
                // only if the counts gain in every bucket is at least 1 (if not zero)
                double sig2ContribGain = 0;

                for (int b = 0; b < mBucketCount; ++b)
                {
                    // if(b == eb) // looking to switch a different bucket from the exhausted one?
                    //    continue;

                    if (cbData[b][s2] == 0)
                        continue;

                    // the first sig can be reduced to zero, but the second sig has a limit to how much it can gain
                    // want to find the maximum that the
                    double countExcess = mCounts[b] - lessSig1Counts[b];

                    if (countExcess < 1)
                    {
                        // reducing sig's contribution in the exhausted bucket had no effect on this sig
                        sig2ContribGain = 0;
                        break;
                    }

                    double alloc = min(countExcess * cbData[b][s2], mCountsTotal);

                    if (sig2ContribGain == 0 || alloc < sig2ContribGain)
                    {
                        sig2ContribGain = alloc;
                    }
                }

                if (sig2ContribGain == 0)
                    continue;

                // check that this gain translates to at least a 1-lot increase in every bucket
                for (int b = 0; b < mBucketCount; ++b)
                {
                    // if(b == eb) // looking to switch a different bucket from the exhausted one?
                    //    continue;

                    if (sigData[b][s2] == 0)
                        continue;

                    if(lessThan(sig2ContribGain * sigData[b][s2], 1))
                    {
                        sig2ContribGain = 0;
                        break;
                    }
                }

                if (sig2ContribGain > topSig2ContribGain)
                {
                    topSig2 = s2;
                    topSig2ContribGain = sig2ContribGain;
                }
            }

            if(topSig2ContribGain > cumulativeSig1ContribLoss)
            {
                // validate the proposed change
                if(greaterThan(mContribTotal + topSig2ContribGain - cumulativeSig1ContribLoss, mCountsTotal))
                {
                    double maxIncrease = mCountsTotal - (mContribTotal - cumulativeSig1ContribLoss);
                    topSig2ContribGain = min(topSig2ContribGain, maxIncrease);
                }

                if(mLogVerbose)
                {
                    LOGGER.debug(String.format("reduceSig(%s cur=%s loss=%s) gainSig(%d cur=%s gain=%s)",
                            rs, sizeToStr(mContribs[rs]), sizeToStr(cumulativeSig1ContribLoss),
                            topSig2, sizeToStr(mContribs[topSig2]), sizeToStr(topSig2ContribGain)));
                }

                applyContribution(rs, -cumulativeSig1ContribLoss);
                applyContribution(topSig2, topSig2ContribGain);
                logStats();

                totalSig1ContribLoss += cumulativeSig1ContribLoss;
                cumulativeSig1ContribLoss = 0;
            }
            else
            {
                // otherwise reduce sig 1's contribution again and keep searching
            }
        }

        return totalSig1ContribLoss > 0;
    }

    private void calcAllContributions()
    {
        for(int s = 0; s < mSigCount; ++s)
        {
            if(mZeroedSigs.contains(s))
                continue;

            calcSigContribution(s);
        }
    }

    private void calcSigContribution(int sig)
    {
        double[][] cbData = mSigCostBasis.getData();

        double minAlloc = 0;

        for (int b = 0; b < mBucketCount; ++b)
        {
            if (cbData[b][sig] == 0)
                continue;

            if (mCurrentCounts[b] >= mCounts[b])
            {
                minAlloc = 0;
                break;
            }

            double alloc = (mCounts[b] - mCurrentCounts[b]) * cbData[b][sig];

            alloc = min(alloc, mCountsTotal);

            if (minAlloc == 0 || alloc < minAlloc)
            {
                minAlloc = alloc;
            }
        }

        if(minAlloc == 0)
            return;

        applyContribution(sig, minAlloc);
    }

    private void applyContribution(int sig, double newContrib)
    {
        // check exceeding total or going negative
        if(greaterThan(mContribTotal + newContrib, mCountsTotal) || lessThan(mContribTotal + newContrib, 0))
        {
            LOGGER.error(String.format("invalid contrib addition(%.1f + %.1f)", mContribTotal, newContrib));
            mIsValid = false;
            return;
        }

        if(mContribs[sig] + newContrib < 0)
        {
            LOGGER.error(String.format("invalid sig(%d) contrib reduction(%.1f + %.1f)", sig, mContribs[sig], newContrib));
            mIsValid = false;
            return;
        }

        mContribs[sig] += newContrib;
        mContribTotal += newContrib;

        final double[][] sigData = mSigs.getData();

        for (int b = 0; b < mBucketCount; ++b)
        {
            double newCount = newContrib * sigData[b][sig];

            if(greaterThan(mCurrentCounts[b] + newCount, mCounts[b]))
            {
                mIsValid = false;
                return;
            }

            mCurrentCounts[b] += newCount;
        }
    }

    private void zeroSigContrib(int sig)
    {
        if(mZeroedSigs.contains(sig))
        {
            LOGGER.error("sample({}) sig({}) previously zeroed", mSampleId, sig);
            mIsValid = false;
            return;
        }

        mZeroedSigs.add(sig);

        double[][] sigData = mSigs.getData();
        double[][] cbData = mSigCostBasis.getData();

        for (int b = 0; b < mBucketCount; ++b)
        {
            mCurrentCounts[b] -= mContribs[sig] * sigData[b][sig];

            if(lessThan(mCurrentCounts[b], 0))
            {
                LOGGER.error(String.format("invalid sig(%d) contrib(%.1f) zero reduction for currentCount(%.1f) and sigRatio(%.4f)",
                        sig, mContribs[sig], mCurrentCounts[b], sigData[b][sig]));

                mIsValid = false;
                return;
            }
        }

        mContribTotal -= mContribs[sig];
        mContribs[sig] = 0;

        for (int b = 0; b < mBucketCount; ++b)
        {
            sigData[b][sig] = 0;
            cbData[b][sig] = 0;
        }
    }

    private void recalcStats()
    {
        // validate current state
        calcResiduals();

        mContribTotal = sumVector(mContribs);

        for(int i = 0; i < mBucketCount; ++i)
        {
            if(greaterThan(mCurrentCounts[i], mCounts[i]))
            {
                LOGGER.error(String.format("sample(%d) bucket currentCount(%.1f) exceeds count(%.1f)", mSampleId, mCurrentCounts[i], mCounts[i]));
                mIsValid = false;
            }
        }

        if(greaterThan(mContribTotal, mCountsTotal))
        {
            LOGGER.error(String.format("sample(%d) contribTotal(%.1f) exceeds totalCount(%.1f)", mSampleId, mContribTotal, mCountsTotal));
            mIsValid = false;
        }

        if(abs(mCountsTotal - mContribTotal - mResiduals) >= 1)
        {
            LOGGER.error(String.format("sample(%d) totalCount(%.1f) less contribTotal(%.1f) != residuals(%.1f))", mSampleId, mCountsTotal, mContribTotal, mResiduals));
            mIsValid = false;
        }

        List<Integer> sortedContribIndices = getSortedVectorIndices(mContribs, false);

        mHasLowContribSigs = false;

        double allocAboveMinRequired = 0;

        for (Integer s : sortedContribIndices)
        {
            if (mContribs[s] == 0)
                continue;

            if(mContribs[s] / mVarCount >= mMinContribPercent)
                allocAboveMinRequired += mContribs[s];

            if(s == mRequiredSig)
                continue;

            if(mContribs[s] / mVarCount < mMinContribPercent)
                mHasLowContribSigs = true;
        }

        mCurrentAllocPerc = min(allocAboveMinRequired/mVarCount, 1);
        mIsFullyAllocated = mCurrentAllocPerc >= mTargetAllocPercent;
    }

    private void clearLowContribSigs()
    {
        if(!mIsFullyAllocated)
            return;

        for (int s = 0; s < mSigCount; ++s)
        {
            if(mZeroedSigs.contains(s))
                continue;

            if (mContribs[s] / mVarCount < mMinContribPercent)
            {
                zeroSigContrib(s);
            }
        }

        recalcStats();
    }

    private void logStats()
    {
        recalcStats();

        if(!mLogVerbose)
            return;

        double actualResiduals = max(mVarCount - mContribTotal, 0);

        LOGGER.debug(String.format("sample(%d) totalCount(%s wn=%s) contrib(%s perc=%.3f) residuals(%s perc=%.3f)",
                mSampleId, sizeToStr(mVarCount), sizeToStr(mCountsTotal), sizeToStr(mContribTotal), mCurrentAllocPerc,
                sizeToStr(actualResiduals), actualResiduals/mVarCount));

        String contribStr = "";
        int sigCount = 0;

        List<Integer> sortedContribIndices = getSortedVectorIndices(mContribs, false);

        for (Integer s : sortedContribIndices)
        {
            if(mContribs[s] == 0)
                break;

            ++sigCount;

            if(!contribStr.isEmpty())
            {
                contribStr += ", ";
            }

            contribStr += String.format("%d = %s perc=%.3f", s, sizeToStr(mContribs[s]), min(mContribs[s]/mVarCount, 1));
        }

        if(!contribStr.isEmpty())
            LOGGER.debug("sample({}) sigs({}) contribs: {}", mSampleId, sigCount, contribStr);
    }

    private List<Integer> getExhaustedBuckets()
    {
        List<Integer> exhaustedBuckets = Lists.newArrayList();

        for (int b = 0; b < mBucketCount; ++b)
        {
            if(!lessThan(mCurrentCounts[b], mCounts[b]))
                exhaustedBuckets.add(b);
        }

        return exhaustedBuckets;
    }

    private void calcResiduals()
    {
        mResiduals = 0;

        for(int i = 0; i < mBucketCount; ++i)
        {
            mResiduals += abs(mCounts[i] - mCurrentCounts[i]);
        }
    }
}
