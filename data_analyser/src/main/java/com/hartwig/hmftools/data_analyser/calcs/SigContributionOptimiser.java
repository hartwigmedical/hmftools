package com.hartwig.hmftools.data_analyser.calcs;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.capValue;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.copyVector;
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

    private final double[] mCounts;
    private final double[] mCountsMargin;
    private final List<double[]> mRatiosCollection;

    private final double[] mInitContribs;

    private double[] mCurrentCounts;
    private double[] mContribs;
    private double[] mMarginCounts;

    private NmfMatrix mSigs;
    private NmfMatrix mSigCostBasis;

    private int mBucketCount;
    private int mSigCount;
    private double mTotalCount;
    private double mResiduals;
    private double mCurrentAllocPerc = 0;
    private double mInitAllocPerc = 0;
    private double mMinContribChange;

    // constants and control config
    private static double MIN_COUNT_CHG_PERC = 0.001;
    private static double MIN_COUNT_CHG = 1;

    private static final Logger LOGGER = LogManager.getLogger(SigContributionOptimiser.class);

    public SigContributionOptimiser(
            int sampleId, final double[] counts, final double[] countsMargin, final List<double[]> ratiosCollection,
            final double[] contribs)
    {
        mSampleId = sampleId;

        mCounts = counts;
        mCountsMargin = countsMargin;

        mRatiosCollection = Lists.newArrayList();

        mSigCount = ratiosCollection.size();
        mBucketCount = counts.length;

        mCurrentCounts = new double[mBucketCount];
        mMarginCounts = new double[mBucketCount];
        copyVector(mCounts, mMarginCounts);
        sumVectors(mCountsMargin, mMarginCounts);

        mRatiosCollection.addAll(ratiosCollection);

        mContribs = new double[mSigCount];
        mInitContribs = new double[mSigCount];
        copyVector(contribs, mContribs);
        copyVector(contribs, mInitContribs);

        mTotalCount = sumVector(mCounts);
        mResiduals = 0;
        mCurrentAllocPerc = 0;
        mInitAllocPerc = 0;

        mMinContribChange = max(mTotalCount * MIN_COUNT_CHG_PERC, MIN_COUNT_CHG);

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
    public double getResiduals() { return mResiduals; }
    public double getAllocPerc() { return mCurrentAllocPerc; }
    public boolean isValid() { return mIsValid; }

    private void applyInitialContributions()
    {
        double[][] sigData = mSigs.getData();
    
        for(int s = 0; s < mSigCount; ++s)
        {
            double minAlloc = 0;

            for (int b = 0; b < mBucketCount; ++b)
            {
                if (sigData[b][s] == 0)
                    continue;

                if (mCurrentCounts[b] >= mCounts[b])
                {
                    minAlloc = 0;
                    break;
                }

                double alloc = (mCounts[b] - mCurrentCounts[b]) / sigData[b][s];

                if (minAlloc == 0 || alloc < minAlloc)
                {
                    minAlloc = alloc;
                }
            }

            if(minAlloc == 0)
                continue;

            applyContribution(s, minAlloc);
        }

        double initContrib = sumVector(mInitContribs);
        mInitAllocPerc = initContrib/mTotalCount;
    }

    public boolean fitToSample(double targetAllocPercent)
    {
        // first apply all the sigs' allocations in turn to exhaust as much of the mCounts as possible
        applyInitialContributions();
        logStats();

        // continue optimising until there has been a run through all exhausted buckets with no improvements
        boolean fullyAllocated = mCurrentAllocPerc >= targetAllocPercent;

        int iterations = 0;
        int maxIterations = 50;

        while(iterations < maxIterations)
        {
            List<Integer> exhaustedBuckets = getExhaustedBuckets();
            List<Integer> sortedContribIndices = getSortedVectorIndices(mContribs, false);

            boolean foundTransfers = false;

            for (Integer eb : exhaustedBuckets)
            {
                for (Integer sig : sortedContribIndices)
                {
                    foundTransfers = reduceSigContribution(sig, eb);

                    if (foundTransfers)
                    {
                        // if an contributions have changed, need to break and begin the process again
                        logStats();
                        break;
                    }
                }

                if (foundTransfers)
                    break;
            }

            if(mCurrentAllocPerc >= targetAllocPercent)
            {
                LOGGER.debug(String.format("sample(%d) allocPerc(%.3f) exceeds target", mSampleId, mCurrentAllocPerc));
                break;
            }

            ++iterations;
        }

        return true;
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
        // double[] ebSigCounts = new double[mSigCount];
        // double[] ebSigCostBasis = new double[mSigCount];

        double[][] sigData = mSigs.getData();
        double[][] cbData = mSigCostBasis.getData();

        double[] lessSig1Counts = new double[mBucketCount];

        // calc the minimum the contribution loss in the exhausted bucket for the sig being reduced
        double sig1EBContribLoss = mMinContribChange / cbData[eb][rs];

        if (sig1EBContribLoss == 0)
            return false;

        double cumulativeSig1ContribLoss = 0;
        double totalSig1ContribLoss = 0;

        while (lessSig1Counts[eb] >= 0)
        {
            // cache the best option for transfer
            int topSig2 = 0;
            double topSig2ContribGain = 0;

            copyVector(mCurrentCounts, lessSig1Counts);

            cumulativeSig1ContribLoss += sig1EBContribLoss;

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

                    if (sigData[b][s2] == 0)
                        continue;

                    // the first sig can be reduced to zero, but the second sig has a limit to how much it can gain
                    // want to find the maximum that the
                    double countExcess = mCounts[b] - lessSig1Counts[b];

                    if (countExcess < mMinContribChange)
                    {
                        // reducing sig's contribution in the exhausted bucket had no effect on this sig
                        sig2ContribGain = 0;
                        break;
                    }

                    double alloc = min(countExcess / sigData[b][s2], mTotalCount);

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

                    if(sig2ContribGain * sigData[b][s2] < mMinContribChange)
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
                LOGGER.debug(String.format("reduceSig(%s loss=%s) gainSig(%d gain=%s)",
                        rs, sizeToStr(cumulativeSig1ContribLoss), topSig2, sizeToStr(topSig2ContribGain)));

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

    private void applyContribution(int sig, double newContrib)
    {
        if((mContribs[sig] + newContrib > mTotalCount) || (mContribs[sig] + newContrib < 0)) // check exceeding total or going negative
        {
            LOGGER.error(String.format("invalid contrib adjustment(%.1f + %.1f)", mContribs[sig], newContrib));
        }

        mContribs[sig] += newContrib;

        for (int b = 0; b < mBucketCount; ++b)
        {
            double newCount = newContrib * mSigs.get(b, sig);
            mCurrentCounts[b] = min(mCurrentCounts[b] + newCount, mCounts[b]);
        }
    }


    private void logStats()
    {
        calcResiduals();
        double currentContrib = sumVector(mContribs);
        mCurrentAllocPerc = currentContrib/mTotalCount;

        LOGGER.debug(String.format("sample(%d) totalCount(%s) contrib(%s perc=%.3f init=%.3f) residuals(%s perc=%.3f)",
                mSampleId, sizeToStr(mTotalCount), sizeToStr(currentContrib), mCurrentAllocPerc, mInitAllocPerc,
                sizeToStr(mResiduals), mResiduals/mTotalCount));

        String contribStr = "";
        for(int s = 0; s < mSigCount; ++s)
        {
            if(s > 0)
            {
                contribStr += ", ";
            }

            contribStr += String.format("%d = %s perc=%.3f", s, sizeToStr(mContribs[s]), mContribs[s]/mTotalCount);
        }

        LOGGER.debug("sigs({}) contribs: {}", mSigCount, contribStr);
    }

    private List<Integer> getExhaustedBuckets()
    {
        List<Integer> exhausedBuckets = Lists.newArrayList();

        for (int b = 0; b < mBucketCount; ++b)
        {
            if(mCurrentCounts[b] >= mCounts[b])
                exhausedBuckets.add(b);
        }

        return exhausedBuckets;
    }

    private void calcResiduals()
    {
        mResiduals = 0;

        for(int i = 0; i < mBucketCount; ++i)
        {
            mResiduals += abs(mCurrentCounts[i] - mCounts[i]);
        }
    }


}
