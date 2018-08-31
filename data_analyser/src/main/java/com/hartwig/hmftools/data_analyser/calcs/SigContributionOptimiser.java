package com.hartwig.hmftools.data_analyser.calcs;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.calcCSS;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.copyVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.doubleToStr;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.greaterThan;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.initVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.lessThan;
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
    private int[] mSigIds;

    private double[] mInitContribs;
    private double mMinContribPercent; // each sig's final contrib must be above this level as percent of total variants
    private int mMinContribCount; // each sig's final contrib must be above this level in absolute terms
    private double mTargetAllocPercent; // target total allocation to exit the fit
    private int mTargetSig; // to check if a specific sig remains above the require contribution percent
    private int mRequiredSig; // keep this specific sig even if it falls below the require contribution percent
    private List<Integer> mZeroedSigs; // the signature bucket ratios

    private double[] mCounts; // sample counts, optionally including noise
    private double[] mCurrentCounts;
    private double[] mContribs;
    private double[] mMaxContribs; // per sig if by itself
    private double mContribTotal;

    private NmfMatrix mSigs;
    private NmfMatrix mSigCostBasis;

    private int mBucketCount;
    private int mSigCount;
    private double mCountsTotal; // total of the counts including any applied noise
    private double mVarCount; // total of the actual variant counts
    private double mResiduals;
    private double mCurrentAllocTotal; // sum of contributions above the require percent and capped at actual counts, not noise
    private double mCurrentAllocPerc;
    private double mInitAllocPerc;
    private double mMinContribChange;
    private boolean mHasLowContribSigs;
    private boolean mIsFullyAllocated;

    private int mNoImproveCount;
    private double mLastImprovedAllocPercent;

    private boolean mLogVerbose;
    private boolean mApplyNoise;

    // constants and control config
    private static double MIN_COUNT_CHG_PERC = 0.001;
    private static double MIN_COUNT_CHG = 1;

    private static final Logger LOGGER = LogManager.getLogger(SigContributionOptimiser.class);

    public SigContributionOptimiser(int bucketCount, boolean logVerbose, double targetAllocPercent, boolean applyNoise)
    {
        mBucketCount = bucketCount;

        mLogVerbose = logVerbose;
        mMinContribPercent = 0.001;
        mTargetAllocPercent = targetAllocPercent;
        mApplyNoise = applyNoise;

        mTargetSig = -1;
        mRequiredSig = -1;
        mZeroedSigs = Lists.newArrayList();
        mRatiosCollection = Lists.newArrayList();

        mRawCounts = new double[mBucketCount];
        mCurrentCounts = new double[mBucketCount];
        mCounts = new double[mBucketCount];
        mCountsNoise = new double[mBucketCount];
    }

    public void initialise(int sampleId, final double[] counts, final double[] countsNoise, final List<double[]> ratiosCollection, double minSigPercent, int minAllocCount)
    {
        mSampleId = sampleId;
        mMinContribPercent = minSigPercent;
        mMinContribCount = minAllocCount;

        mSigCount = ratiosCollection.size();
        mSigIds = new int[mSigCount];
        for(int i = 0; i < mSigCount; ++i)
        {
            mSigIds[i] = i;
        }

        copyVector(counts, mRawCounts);
        copyVector(countsNoise, mCountsNoise);
        copyVector(counts, mCounts);

        mTargetSig = -1;
        mRequiredSig = -1;
        mZeroedSigs.clear();

        sumVectors(mCountsNoise, mCounts);

        mRatiosCollection.clear();
        mRatiosCollection.addAll(ratiosCollection);

        mContribs = new double[mSigCount];
        mInitContribs = new double[mSigCount];
        mMaxContribs = new double[mSigCount];
        mContribTotal = 0;

        mVarCount = sumVector(mRawCounts);
        mCountsTotal = sumVector(mCounts);
        mResiduals = 0;
        mCurrentAllocPerc = 0;
        mInitAllocPerc = 0;
        mHasLowContribSigs = false;
        mIsFullyAllocated = false;
        mNoImproveCount = 0;
        mLastImprovedAllocPercent = 0;

        calcMinContribChange();

        for (int b = 0; b < mBucketCount; ++b)
        {
            mCurrentCounts[b] = 0;
        }

        // extract sigs their cost basis (just the inverse)
        mSigs = new NmfMatrix(mBucketCount, mSigCount);
        mSigCostBasis = new NmfMatrix(mBucketCount, mSigCount);
        // mSigCss = new NmfMatrix(mSigCount, mSigCount);

        if (mContribs.length != ratiosCollection.size())
        {
            mIsValid = false;
            return;
        }

        double[][] sigData = mSigs.getData();
        // double[][] cssData = mSigCss.getData();
        double[][] cbData = mSigCostBasis.getData();

        for (int sig = 0; sig < mSigCount; ++sig)
        {
            double[] sigRatios = ratiosCollection.get(sig);

            if (sigRatios.length != mBucketCount)
            {
                mIsValid = false;
                return;
            }

            for (int b = 0; b < mBucketCount; ++b)
            {
                sigData[b][sig] = sigRatios[b];

                if (sigRatios[b] > 0)
                    cbData[b][sig] = 1 / sigRatios[b];
            }

            mContribs[sig] = 0;
        }

        mIsValid = true;
    }

    public void setSigIds(List<Integer> sigIds)
    {
        if(sigIds.size() != mSigCount)
            return;

        for(int sig = 0; sig < mSigCount; ++sig)
        {
            mSigIds[sig] = sigIds.get(sig);
        }
    }

    public final double[] getFittedCounts()
    {
        return mCurrentCounts;
    }

    public final double[] getContribs()
    {
        return mContribs;
    }

    public final double getContribTotal()
    {
        return mContribTotal;
    }

    public double getAllocPerc()
    {
        return mCurrentAllocPerc;
    }

    public boolean isValid()
    {
        return mIsValid;
    }

    public boolean isFullyAllocated()
    {
        return mIsFullyAllocated;
    }

    public void setTargetSig(int sig)
    {
        mTargetSig = sig;
    }

    public void setRequiredSig(int sig)
    {
        mRequiredSig = sig;
    }

    public void setLogVerbose(boolean toggle)
    {
        mLogVerbose = toggle;
    }

    private void applyInitialContributions()
    {
        // calc the theoretical max for each sig
        /*
        for (int s = 0; s < mSigCount; ++s)
        {
            mMaxContribs[s] = calcSigContribution(s);
        }
        */

        // try various approaches to find the best init allocations
        double[] sigContribGains = new double[mSigCount];
        List<Integer> emptyBuckets = Lists.newArrayList();
        double initSigLoss = testSigReduction(-1, emptyBuckets, sigContribGains);
        double initGain = sumVector(sigContribGains);

        if(initGain > 0)
        {
            applySigAdjustments(-1, 0, sigContribGains);
        }

        copyVector(mContribs, mInitContribs);
        double initContrib = sumVector(mInitContribs);
        mInitAllocPerc = initContrib / mVarCount;

        for (int s = 0; s < mSigCount; ++s)
        {
            if (s == mRequiredSig)
                continue;

            if (!aboveMinReqContrib(mInitContribs[s]))
                mHasLowContribSigs = true;
        }
    }

    public boolean fitToSample()
    {
        if (!mIsValid)
            return false;

        // first apply all the sigs' allocations in turn to exhaust as much of the mCounts as possible
        applyInitialContributions();
        logStats();

        if (!mIsValid)
            return false;

        if (mIsFullyAllocated)
        {
            clearLowContribSigs();
            return true;
        }

        boolean foundImprovements = findAdjustments();

        while (foundImprovements || mHasLowContribSigs)
        {
            // strip out the worst contributor if below teh required threshold and try again
            List<Integer> sortedContribIndices = getSortedVectorIndices(mContribs, true);

            for (Integer s : sortedContribIndices)
            {
                if (mContribs[s] == 0)
                    continue;

                if (s == mRequiredSig)
                    continue;

                if (!aboveMinReqContrib(mContribs[s]))
                {
                    zeroSigContrib(s);

                    if (s == mTargetSig)
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
            foundImprovements = findAdjustments();

            if (!mIsValid)
                return false;

            if(mNoImproveCount > 10)
            {
                foundImprovements = false;
            }

            if (mIsFullyAllocated)
            {
                clearLowContribSigs();
                return true;
            }
        }

        return mIsValid;
    }

    private boolean findAdjustments()
    {
        // continue optimising until there has been a run through all exhausted buckets with no improvements
        int iterations = 0;
        int maxIterations = 100;

        boolean recalcRequired = false;

        while (iterations < maxIterations)
        {
            List<Integer> exhaustedBuckets = getExhaustedBuckets(0.9);

            if (exhaustedBuckets.isEmpty())
            {
                calcAllContributions();
                exhaustedBuckets = getExhaustedBuckets(0.9);

                if (exhaustedBuckets.isEmpty())
                    return false;
            }

            boolean foundTransfers = false;

            /*
            List<Integer> sortedContribIndices = getSortedVectorIndices(mContribs, false);

            for (Integer sig : sortedContribIndices)
            {
                foundTransfers = reduceSigContribution(sig, exhaustedBuckets);

                if(!mIsValid)
                    return false;

                if (foundTransfers)
                {
                    // if any contributions have changed, need to break and begin the process again
                    recalcRequired = true;
                    logStats();
                    break;
                }
            }
            */

            double[] maxOtherSigContribGains = new double[mSigCount];
            double maxReducedSigLoss = 0;
            double maxNetGain = 0;
            int maxReducedSig = -1;

            // find the sig with the max net gain from being reduced
            for (int sig = 0; sig < mSigCount; ++sig)
            {
                double[] otherSigContribGains = new double[mSigCount];
                double reducedSigContribLoss = testSigReduction(sig, exhaustedBuckets, otherSigContribGains);
                double netGain = sumVector(otherSigContribGains) - reducedSigContribLoss;

                if (reducedSigContribLoss > 0 && netGain > maxNetGain)
                {
                    maxNetGain = netGain;
                    maxReducedSig = sig;
                    maxReducedSigLoss = reducedSigContribLoss;
                    copyVector(otherSigContribGains, maxOtherSigContribGains);
                }
            }

            if (maxReducedSig >= 0)
            {
                applySigAdjustments(maxReducedSig, maxReducedSigLoss, maxOtherSigContribGains);

                if (!mIsValid)
                    return false;

                recalcRequired = true;
                foundTransfers = true;
                logStats();
                break;
            }

            if (mIsFullyAllocated)
            {
                if (mLogVerbose)
                    LOGGER.debug(String.format("sample(%d) allocPerc(%.3f) exceeds target", mSampleId, mCurrentAllocPerc));

                break;
            }

            if (!foundTransfers) // no point continuing to look
                break;

            ++iterations;
        }

        return recalcRequired;
    }

    private double testSigReduction(int rs, List<Integer> exhaustedBuckets, double[] otherSigContribGains)
    {
        double[][] sigData = mSigs.getData();
        double[][] cbData = mSigCostBasis.getData();

        double sig1ContribLoss = 0;
        double[] maxOtherSigContribs = new double[mSigCount];

        double[] lessSig1Counts = new double[mBucketCount];
        copyVector(mCurrentCounts, lessSig1Counts);

        boolean initialTest = exhaustedBuckets.isEmpty();

        if(!initialTest)
        {
            // search through the exhausted buckets to find the one which delivers the largest gain in the reducing sig
            double minSig1ContribLoss = 0;
            for (Integer eb : exhaustedBuckets)
            {
                if (cbData[eb][rs] == 0)
                    continue; // return 0;

                // calc the minimum the contribution loss in the exhausted bucket for the sig being reduced
                double contribLoss = mMinContribChange * cbData[eb][rs];

                if (minSig1ContribLoss == 0 || contribLoss < minSig1ContribLoss)
                {
                    minSig1ContribLoss = contribLoss;
                }
            }

            if (minSig1ContribLoss == 0)
                return 0;

            sig1ContribLoss = min(minSig1ContribLoss, mContribs[rs]); // cannot reduce past zero

            // remove this sig across the board from a 1-lot contrib to this exhausted bucket
            for (int b = 0; b < mBucketCount; ++b)
            {
                lessSig1Counts[b] -= sig1ContribLoss * sigData[b][rs];

                if (lessThan(lessSig1Counts[b], 0))
                {
                    mIsValid = false;
                    LOGGER.error(String.format("bucket(%d) currentCount(%.1f) reduced below zero from sig(%d contrib=%.1f)", b, lessSig1Counts[b], mSigIds[rs], sig1ContribLoss));
                    return 0;
                }
            }
        }

        // now look at the potential gain to all the sigs in each bucket, having had sig1 removed
        for (int s2 = 0; s2 < mSigCount; ++s2)
        {
            if (rs == s2)
                continue;

            double minAlloc = calcSigContribution(s2, lessSig1Counts);

            if (minAlloc > 0)
            {
                maxOtherSigContribs[s2] = minAlloc;
            }
        }

        double totalOtherSigsGain = sumVector(maxOtherSigContribs);

        if (totalOtherSigsGain < sig1ContribLoss * 0.5)
            return 0;

        // sort into order of greatest effect
        // List<Integer> otherSigContribIndices = getSortedVectorIndices(maxOtherSigGains, false);

        // test out the proposed change to find the max that can be applied
        double[] testSig1Counts = new double[mBucketCount];
        double[] testOtherSigContribs = new double[mSigCount];
        double maxOtherSigsGain = 0;

        // runs 0 and 1, sort ascending then descending with total allocation top-down each time
        // runs 2 and 3, sort ascending then descending with partial allocation top-down each time
        // runs 4 and 5, sort descending with BG first, with full then partial allocation
        for (int i = 0; i < 6; ++i)
        {
            copyVector(lessSig1Counts, testSig1Counts);

            // sort descending except on first iteration
            List<Integer> otherSigContribIndices = getSortedVectorIndices(maxOtherSigContribs, (i == 0 || i == 2));
            // List<Integer> otherSigContribIndices = getSortedVectorIndices(maxOtherSigContribs, false);

            // try the likely order from the original discovery fit process
            if(i >= 4)
            {
                if(!initialTest)
                    break;

                if(otherSigContribIndices.get(0) == mRequiredSig)
                    break;

                // put the required sig first
                otherSigContribIndices.add(0, mRequiredSig);
            }

            if (i > 0)
                initVector(testOtherSigContribs, 0);

            double allocPerc = (i == 2 || i == 3 || i == 5) ? 0.75 : 1;
            boolean foundAdjusts = true;
            int iterations = 0;
            int maxIteration = 5;

            while (foundAdjusts && iterations < maxIteration)
            {
                foundAdjusts = false;

                for (Integer s2 : otherSigContribIndices)
                {
                    if (s2 == rs || maxOtherSigContribs[s2] == 0)
                        continue;

                    // LOGGER.debug(String.format("testing sig(%d) gain(%s) vs current(%s)", s2, sizeToStr(otherSigContribGains[s2]), sizeToStr(mContribs[s2])));

                    double newAlloc = calcSigContribution(s2, testSig1Counts);

                    if (newAlloc <= 0)
                        continue;

                    if (iterations < maxIteration - 1)
                        newAlloc *= allocPerc;

                    testOtherSigContribs[s2] += newAlloc; // take adjustment if required
                    foundAdjusts = true;

                    for (int b = 0; b < mBucketCount; ++b)
                    {
                        double newCount = newAlloc * sigData[b][s2];

                        if(newCount == 0)
                            continue;

                        if (greaterThan(testSig1Counts[b] + newCount, mCounts[b]))
                        {
                            testOtherSigContribs[s2] = 0; // suggests an error in the calcSigContrib function
                            break;
                        }

                        testSig1Counts[b] += newCount;
                    }
                }

                ++iterations;
            }

            // finally factor in potentially being able to add back in some portion of the sig that was reduced
            if(!initialTest)
            {
                double newAlloc = calcSigContribution(rs, testSig1Counts);

                if (newAlloc > 0)
                {
                    testOtherSigContribs[rs] += newAlloc;
                }
            }

            double sigContribsTotal = sumVector(testOtherSigContribs);

            if (sigContribsTotal > maxOtherSigsGain)
            {
                // take the top allocation combination
                maxOtherSigsGain = sigContribsTotal;
                copyVector(testOtherSigContribs, otherSigContribGains);
            }
        }

        if (maxOtherSigsGain < sig1ContribLoss)
            return 0;

        return sig1ContribLoss;
    }

    private void applySigAdjustments(int rs, double sigContribLoss, double[] otherSigContribGains)
    {
        if(rs >= 0)
        {
            if (mLogVerbose)
            {
                double totalActualOtherGain = sumVector(otherSigContribGains);

                LOGGER.debug(String.format("reduceSig(%s cur=%s loss=%s) for other sigs gain(%s)",
                        mSigIds[rs], doubleToStr(mContribs[rs]), doubleToStr(sigContribLoss), doubleToStr(totalActualOtherGain)));
            }

            applyContribution(rs, -sigContribLoss);
        }

        List<Integer> otherSigContribIndices = getSortedVectorIndices(otherSigContribGains, false);

        for (Integer s2 : otherSigContribIndices)
        {
            if (s2 == rs)
                continue;

            if (otherSigContribGains[s2] == 0)
                break;

            if (mLogVerbose)
            {
                LOGGER.debug(String.format("apply sig(%d) contrib(%s -> %s gain=%s)",
                        mSigIds[s2], doubleToStr(mContribs[s2]), doubleToStr(mContribs[s2] + otherSigContribGains[s2]), doubleToStr(otherSigContribGains[s2])));
            }

            applyContribution(s2, otherSigContribGains[s2]);
        }

        // logStats();
    }

    private void calcAllContributions()
    {
        for (int s = 0; s < mSigCount; ++s)
        {
            if (mZeroedSigs.contains(s))
                continue;

            double alloc = calcSigContribution(s);
            applyContribution(s, alloc);
        }
    }

    private double calcSigContribution(int sig)
    {
        return calcSigContribution(sig, mCurrentCounts);
    }

    private double calcSigContribution(int sig, final double[] currentCounts)
    {
        double[][] cbData = mSigCostBasis.getData();

        double minAlloc = 0;

        for (int b = 0; b < mBucketCount; ++b)
        {
            if (cbData[b][sig] == 0)
                continue;

            if (currentCounts[b] >= mCounts[b])
            {
                minAlloc = 0;
                break;
            }

            double alloc = (mCounts[b] - currentCounts[b]) * cbData[b][sig];

            alloc = min(alloc, mCountsTotal);

            if (minAlloc == 0 || alloc < minAlloc)
            {
                minAlloc = alloc;
            }
        }

        return minAlloc;
    }

    private void applyContribution(int sig, double newContrib)
    {
        if(newContrib == 0)
            return;

        // check exceeding total or going negative
        if (greaterThan(mContribTotal + newContrib, mCountsTotal) || lessThan(mContribTotal + newContrib, 0))
        {
            LOGGER.error(String.format("invalid contrib addition(%.1f + %.1f)", mContribTotal, newContrib));
            mIsValid = false;
            return;
        }

        if (lessThan(mContribs[sig] + newContrib, 0))
        {
            LOGGER.error(String.format("invalid sig(%d) contrib reduction(%.1f + %.1f)", mSigIds[sig], mContribs[sig], newContrib));
            mIsValid = false;
            return;
        }

        mContribs[sig] += newContrib;
        mContribTotal += newContrib;

        final double[][] sigData = mSigs.getData();

        for (int b = 0; b < mBucketCount; ++b)
        {
            double newCount = newContrib * sigData[b][sig];

            if (greaterThan(mCurrentCounts[b] + newCount, mCounts[b]))
            {
                LOGGER.error(String.format("sig(%d contrib=%f) count(%f) + newCount(%f) exceeds maxCount(%f) diff(%f)",
                        mSigIds[sig], mContribs[sig], mCurrentCounts[b], newCount, mCounts[b],
                        mCurrentCounts[b] + newCount - mCounts[b]));
                mIsValid = false;
                return;
            }

            mCurrentCounts[b] += newCount;
        }
    }

    private void zeroSigContrib(int sig)
    {
        if (mZeroedSigs.contains(sig))
        {
            LOGGER.error("sample({}) sig({}) previously zeroed", mSampleId, mSigIds[sig]);
            mIsValid = false;
            return;
        }

        mZeroedSigs.add(sig);

        double[][] sigData = mSigs.getData();
        double[][] cbData = mSigCostBasis.getData();

        for (int b = 0; b < mBucketCount; ++b)
        {
            mCurrentCounts[b] -= mContribs[sig] * sigData[b][sig];

            if (lessThan(mCurrentCounts[b], 0))
            {
                LOGGER.error(String.format("invalid sig(%d) contrib(%.1f) zero reduction for currentCount(%.1f) and sigRatio(%.4f)",
                        mSigIds[sig], mContribs[sig], mCurrentCounts[b], sigData[b][sig]));

                mIsValid = false;
                return;
            }
        }

        if (mLogVerbose)
            LOGGER.debug(String.format("sample(%d) sig(%d) contrib(%s) zeroed", mSampleId, mSigIds[sig], doubleToStr(mContribs[sig])));

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

        for (int i = 0; i < mBucketCount; ++i)
        {
            if (greaterThan(mCurrentCounts[i], mCounts[i]))
            {
                LOGGER.error(String.format("sample(%d) bucket currentCount(%.1f) exceeds count(%.1f)", mSampleId, mCurrentCounts[i], mCounts[i]));
                mIsValid = false;
            }
        }

        if (greaterThan(mContribTotal, mCountsTotal))
        {
            LOGGER.error(String.format("sample(%d) contribTotal(%.1f) exceeds totalCount(%.1f)", mSampleId, mContribTotal, mCountsTotal));
            mIsValid = false;
        }

        if (abs(mCountsTotal - mContribTotal - mResiduals) >= 1)
        {
            LOGGER.error(String.format("sample(%d) totalCount(%.1f) less contribTotal(%.1f) != residuals(%.1f))", mSampleId, mCountsTotal, mContribTotal, mResiduals));
            mIsValid = false;
        }

        List<Integer> sortedContribIndices = getSortedVectorIndices(mContribs, false);

        mHasLowContribSigs = false;

        double[] currentCountsNoNoise = new double[mBucketCount];

        final double[][] sigData = mSigs.getData();
        for (Integer s : sortedContribIndices)
        {
            if (mContribs[s] == 0)
                continue;

            if (aboveMinReqContrib(mContribs[s]))
            {
                // also calc an allocation limited by the actual counts, not factoring in noise
                for (int i = 0; i < mBucketCount; ++i)
                {
                    currentCountsNoNoise[i] = min(currentCountsNoNoise[i] + (mContribs[s] * sigData[i][s]), mRawCounts[i]);
                }
            }

            if (s == mRequiredSig)
                continue;

            if (!aboveMinReqContrib(mContribs[s]))
                mHasLowContribSigs = true;
        }

        mCurrentAllocTotal = sumVector(currentCountsNoNoise);

        mCurrentAllocPerc = min(mCurrentAllocTotal / mVarCount, 1);
        mIsFullyAllocated = mCurrentAllocPerc >= mTargetAllocPercent;

        if(mCurrentAllocPerc > mLastImprovedAllocPercent)
        {
            mLastImprovedAllocPercent = mCurrentAllocPerc;
            mNoImproveCount = 0;
        }
        else
        {
            ++mNoImproveCount;
        }

        calcMinContribChange();
    }

    private void calcMinContribChange()
    {
        /*
        // make the min contrib change a function of the current and target alloc
        // to have it move more aggressively when further out and then more fine-grained close to the target
        double targetRemaining = max(mTargetAllocPercent - mCurrentAllocPerc, 0) / mTargetAllocPercent;

        double upperPercent = 0.02;
        double lowerPercent = 0.001;
        double absMinChange = 0.4;
        double requiredChgPerc = lowerPercent + (upperPercent - lowerPercent) * targetRemaining;

        mMinContribChange = max(mVarCount * requiredChgPerc, absMinChange);
        */

        mMinContribChange = max(mVarCount * MIN_COUNT_CHG_PERC, MIN_COUNT_CHG);
    }

    private boolean aboveMinReqContrib(double contrib)
    {
        return (contrib / mVarCount >= mMinContribPercent) && contrib >= mMinContribCount;
    }

    private void clearLowContribSigs()
    {
        if (!mIsFullyAllocated)
            return;

        for (int s = 0; s < mSigCount; ++s)
        {
            if (mZeroedSigs.contains(s))
                continue;

            if (s == mRequiredSig)
                continue;

            if (!aboveMinReqContrib(mContribs[s]))
            {
                zeroSigContrib(s);
            }
        }

        recalcStats();
    }

    private void logStats()
    {
        recalcStats();

        if (!mLogVerbose)
            return;

        double actualResiduals = max(mVarCount - mContribTotal, 0);

        LOGGER.debug(String.format("sample(%d) totalCount(%s wn=%s) allocated(%s act=%.3f can=%.3f) residuals(%s perc=%.3f) minChg(%.1f) hasLow(%s)",
                mSampleId, doubleToStr(mVarCount), doubleToStr(mCountsTotal), doubleToStr(mContribTotal), mContribTotal / mVarCount,
                mCurrentAllocPerc, doubleToStr(actualResiduals), actualResiduals / mVarCount, mMinContribChange, mHasLowContribSigs));

        String contribStr = "";
        int sigCount = 0;

        List<Integer> sortedContribIndices = getSortedVectorIndices(mContribs, false);

        for (Integer s : sortedContribIndices)
        {
            if (mContribs[s] == 0)
                break;

            ++sigCount;

            if (!contribStr.isEmpty())
            {
                contribStr += ", ";
            }

            contribStr += String.format("%d = %s perc=%.3f", mSigIds[s], doubleToStr(mContribs[s]), min(mContribs[s] / mVarCount, 1));
        }

        if (!contribStr.isEmpty())
        {
            String specSigs = "";
            if(mRequiredSig >= 0)
                specSigs = String.format(" req=%d)", mSigIds[mRequiredSig]);
            else if(mTargetSig >= 0)
                specSigs = String.format(" tar=%d)", mSigIds[mTargetSig]);
            else
                specSigs = ")";

            LOGGER.debug("sample({}) sigs({}{} contribs: {}", mSampleId, sigCount, specSigs, contribStr);
        }
    }

    private List<Integer> getExhaustedBuckets(double percFull)
    {
        List<Integer> exhaustedBuckets = Lists.newArrayList();

        for (int b = 0; b < mBucketCount; ++b)
        {
            if (mCurrentCounts[b] >= percFull * mCounts[b])
                exhaustedBuckets.add(b);
        }

        return exhaustedBuckets;
    }

    private void calcResiduals()
    {
        mResiduals = 0;

        for (int i = 0; i < mBucketCount; ++i)
        {
            mResiduals += abs(mCounts[i] - mCurrentCounts[i]);
        }
    }

    // old method
    private boolean reduceSigContribution(int rs, List<Integer> exhaustedBuckets)
    {
        double[][] sigData = mSigs.getData();
        double[][] cbData = mSigCostBasis.getData();

        // search through the exhausted buckets to find the one which delivers the largest gain in the reducing sig
        double maxSig1EBContribLoss = 0;
        for (Integer eb : exhaustedBuckets)
        {
            if (cbData[eb][rs] == 0)
                return false;

            // calc the minimum the contribution loss in the exhausted bucket for the sig being reduced
            double contribLoss = mMinContribChange * cbData[eb][rs]; // mMinContribChange

            if (maxSig1EBContribLoss == 0 || contribLoss < maxSig1EBContribLoss)
            {
                maxSig1EBContribLoss = contribLoss;
            }
        }

        double[] lessSig1Counts = new double[mBucketCount];

        double sig1EBContribLoss = min(maxSig1EBContribLoss, mContribs[rs]); // cannot reduce past zero

        double cumulativeSig1ContribLoss = 0;
        double totalSig1ContribLoss = 0;

        double[] otherSigContribGains = new double[mSigCount];

        while (mContribs[rs] > 0)
        {
            initVector(otherSigContribGains, 0);

            copyVector(mCurrentCounts, lessSig1Counts);

            cumulativeSig1ContribLoss += sig1EBContribLoss;

            if (mContribs[rs] - cumulativeSig1ContribLoss < 0) // cannot reduce further
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

                double minAlloc = calcSigContribution(s2, lessSig1Counts);

                if (minAlloc > 0)
                {
                    otherSigContribGains[s2] = minAlloc;
                }
            }

            double totalOtherSigGain = sumVector(otherSigContribGains);

            if (totalOtherSigGain < cumulativeSig1ContribLoss)
                break;

            // sort into order of greatest effect
            List<Integer> otherSigContribIndices = getSortedVectorIndices(otherSigContribGains, false);

            double totalActualOtherGain = 0;
            int otherSigGainCount = 0;

            // test out the proposed change to find the max that can be applied
            for (Integer s2 : otherSigContribIndices)
            {
                if (s2 == rs)
                    continue;

                // LOGGER.debug(String.format("testing sig(%d) gain(%s) vs current(%s)", s2, sizeToStr(otherSigContribGains[s2]), sizeToStr(mContribs[s2])));

                double newAlloc = calcSigContribution(s2, lessSig1Counts);

                if (newAlloc <= 0)
                {
                    otherSigContribGains[s2] = 0;
                    break;
                }

                ++otherSigGainCount;
                otherSigContribGains[s2] = newAlloc; // take adjustment if required
                totalActualOtherGain += newAlloc;

                for (int b = 0; b < mBucketCount; ++b)
                {
                    double newCount = newAlloc * sigData[b][s2];

                    if (greaterThan(lessSig1Counts[b] + newCount, mCounts[b]))
                        break;

                    lessSig1Counts[b] += newCount;
                }
            }

            if (totalActualOtherGain < cumulativeSig1ContribLoss)
                break;

            if (mLogVerbose)
            {
                LOGGER.debug(String.format("reduceSig(%s cur=%s loss=%s) for %d actual other sigs gain(%s)", rs, doubleToStr(mContribs[rs]), doubleToStr(cumulativeSig1ContribLoss), otherSigGainCount, doubleToStr(totalActualOtherGain)));
            }

            applyContribution(rs, -cumulativeSig1ContribLoss);

            for (Integer s2 : otherSigContribIndices)
            {
                if (s2 == rs)
                    continue;

                if (otherSigContribGains[s2] == 0)
                    break;

                if (mLogVerbose)
                {
                    LOGGER.debug(String.format("apply sig(%d) contrib(%s -> %s gain=%s)", s2, doubleToStr(mContribs[s2] + otherSigContribGains[s2]), doubleToStr(mContribs[s2]), doubleToStr(otherSigContribGains[s2])));
                }

                applyContribution(s2, otherSigContribGains[s2]);
            }

            logStats();

            totalSig1ContribLoss += cumulativeSig1ContribLoss;
            cumulativeSig1ContribLoss = 0;
        }

        return totalSig1ContribLoss > 0;
    }

}