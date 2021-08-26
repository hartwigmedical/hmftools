package com.hartwig.hmftools.sigs.buckets;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.sigs.DataUtils.doubleToStr;
import static com.hartwig.hmftools.common.sigs.DataUtils.doublesEqual;
import static com.hartwig.hmftools.common.sigs.DataUtils.greaterThan;
import static com.hartwig.hmftools.common.sigs.DataUtils.lessThan;
import static com.hartwig.hmftools.common.utils.VectorUtils.copyVector;
import static com.hartwig.hmftools.common.utils.VectorUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.common.utils.VectorUtils.initVector;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVectors;
import static com.hartwig.hmftools.common.utils.VectorUtils.vectorMultiply;
import static com.hartwig.hmftools.sigs.buckets.BaConfig.MAX_NOISE_ALLOC_PERCENT;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Matrix;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class CountsSigContribOptimiser
{
    private int mSampleId;
    private boolean mIsValid;

    private final double[] mRawCounts; // actual sample counts
    private final double[] mCountsNoise; // noise around the counts
    private List<double[]> mRatiosCollection; // the signature bucket ratios
    private int[] mSigIds;

    private double[] mInitContribs;
    private double mMinContribPercent; // each sig's final contrib must be above this level as percent of total variants
    private int mMinContribCount; // each sig's final contrib must be above this level in absolute terms
    private double mTargetAllocPercent; // target total allocation to exit the fit
    private double mTargetResidualsPercent;
    private int mTargetSig; // to check if a specific sig remains above the require contribution percent
    private int mRequiredSig; // keep this specific sig even if it falls below the require contribution percent
    private double mRequiredSigMinContrib;
    private List<Integer> mZeroedSigs; // the signature bucket ratios

    private final double[] mCounts; // sample counts, optionally including noise
    private final double[] mCurrentCounts;
    private double[] mContribs;
    private double mContribTotal;

    private Matrix mSigs;

    private final int mBucketCount;
    private int mSigCount;
    private double mCountsTotal; // total of the counts plus maximum permitted noise allocation
    private double mRawCountsTotal; // total of the actual variant counts
    private double mMaxNoiseTotal;
    private double mResiduals; // standard meaning - absolute value of diff between actual and allocated (including excess)
    private double mCurrentAllocTotal; // sum of contributions above the require percent and capped at actual counts, not noise
    private double mCurrentAllocPerc;
    private double mInitAllocPerc;
    private double mMinContribChange;
    private boolean mHasLowContribSigs;
    private boolean mIsFullyAllocated;

    // calc state to avoid repeated memory alloc
    private final double[] mCurrentCountsNoNoise;
    private final double[] mReducedRefSigCounts; // reduction to counts in the sig being tested
    private final double[] mTestRefSigCounts;
    private final double[] mTestCounts;

    // diagnostics and stats
    private int mIterations;
    private int mInstances;
    private double mAvgIterations;
    private double mAvgPercImprove;
    private List<Double> mRecentAllocPercents;
    private boolean mStagnantAllocChange;

    private boolean mLogVerbose;
    private boolean mLogVerboseOverride;

    // constants and control config
    private static double MIN_COUNT_CHG_PERC = 0.001;
    private static double MIN_COUNT_CHG = 1;
    private static int MAX_TEST_ITERATIONS = 100;
    private static int MAX_NO_IMPROVE_COUNT = 10;
    private static double REQUIRED_SIG_MIN_PERCENT = 0.001;

    private static final Logger LOGGER = LogManager.getLogger(CountsSigContribOptimiser.class);

    public CountsSigContribOptimiser(int bucketCount, boolean logVerbose, double targetAllocPercent)
    {
        mBucketCount = bucketCount;

        mLogVerbose = logVerbose;
        mMinContribPercent = 0.001;
        mTargetAllocPercent = targetAllocPercent;
        mTargetResidualsPercent = (1 - targetAllocPercent) * 5;

        mTargetSig = -1;
        mRequiredSig = -1;
        mRequiredSigMinContrib = 0;
        mZeroedSigs = Lists.newArrayList();
        mRatiosCollection = Lists.newArrayList();
        mRecentAllocPercents = Lists.newArrayList();

        mRawCounts = new double[mBucketCount];
        mCurrentCounts = new double[mBucketCount];
        mCounts = new double[mBucketCount];
        mCountsNoise = new double[mBucketCount];

        mCurrentCountsNoNoise = new double[mBucketCount];
        mReducedRefSigCounts = new double[mBucketCount];
        mTestRefSigCounts = new double[mBucketCount];
        mTestCounts = new double[mBucketCount];

        mInstances = 0;
        mAvgIterations = 0;
        mAvgPercImprove = 0;
    }

    public void initialise(int sampleId, final double[] counts, final double[] noiseCounts, final List<double[]> ratiosCollection,
            double minSigPercent, int minAllocCount)
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
        copyVector(noiseCounts, mCountsNoise);
        copyVector(counts, mCounts);
        initVector(mReducedRefSigCounts, 0);
        initVector(mCurrentCounts, 0);
        initVector(mTestRefSigCounts, 0);
        initVector(mTestCounts, 0);

        mTargetSig = -1;
        mRequiredSig = -1;
        mRequiredSigMinContrib = 0;
        mZeroedSigs.clear();

        sumVectors(mCountsNoise, mCounts);

        mRatiosCollection.clear();
        mRatiosCollection.addAll(ratiosCollection);

        mContribs = new double[mSigCount];
        mInitContribs = new double[mSigCount];
        mContribTotal = 0;

        mRawCountsTotal = sumVector(mRawCounts);
        double noiseTotal = sumVector(noiseCounts);
        mMaxNoiseTotal = min(mRawCountsTotal * MAX_NOISE_ALLOC_PERCENT, noiseTotal);
        mCountsTotal = mRawCountsTotal + mMaxNoiseTotal;

        mResiduals = 0;
        mCurrentAllocPerc = 0;
        mCurrentAllocTotal = 0;
        mInitAllocPerc = 0;
        mHasLowContribSigs = false;
        mIsFullyAllocated = false;
        mRecentAllocPercents.clear();
        mStagnantAllocChange = false;
        mIterations = 0;
        mLogVerboseOverride = false;

        calcMinContribChange();

        // extract sigs their cost basis (just the inverse)
        mSigs = new Matrix(mBucketCount, mSigCount);

        if (mContribs.length != ratiosCollection.size())
        {
            mIsValid = false;
            return;
        }

        double[][] sigData = mSigs.getData();

        for (int sig = 0; sig < mSigCount; ++sig)
        {
            double[] sigRatios = ratiosCollection.get(sig);

            if(!doublesEqual(sumVector(sigRatios), 1))
            {
                LOGGER.error(String.format("sig(%d) has invalid ratios: total(%.6f)", sig, sumVector(sigRatios)));
                mIsValid = false;
                return;
            }

            if (sigRatios.length != mBucketCount)
            {
                mIsValid = false;
                return;
            }

            for (int b = 0; b < mBucketCount; ++b)
            {
                sigData[b][sig] = sigRatios[b];
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

    public final double[] getContribs() { return mContribs; }
    public double getAllocPerc() { return mCurrentAllocPerc; }
    public boolean isValid() { return mIsValid; }
    public int getInstances() { return mInstances; }
    public double getAvgIterations() { return mAvgIterations; }
    public double getAvgImprovePerc() { return mAvgPercImprove; }
    public void setTargetSig(int sig) { mTargetSig = sig; }
    public void setRequiredSig(int sig)
    {
        mRequiredSig = sig;
        mRequiredSigMinContrib = mRawCountsTotal * REQUIRED_SIG_MIN_PERCENT;
    }

    public void setLogVerbose(boolean toggle) { mLogVerbose = toggle; }
    public int contributingSigCount()
    {
        return (int)Arrays.stream(mContribs).filter(x -> x > 0).count();
    }

    public boolean fitToSample()
    {
        if (!mIsValid)
            return false;

        if (!mIsValid)
            return false;

        if (mIsFullyAllocated)
        {
            clearLowContribSigs();
            updateStats();
            return true;
        }

        boolean foundImprovements = findAdjustments();
        logStats();
        ++mIterations;

        boolean targetSigZeroed = false;

        while (foundImprovements || mHasLowContribSigs)
        {
            // strip out the worst contributor if below the required threshold and try again
            List<Integer> sortedContribIndices = getSortedVectorIndices(mContribs, true);

            boolean sigZeroed = false;
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
                        targetSigZeroed = true;
                        break;
                    }

                    //if(mLogVerbose)
                    //    LOGGER.debug("sample({}) removed low-percent sig({})", mSampleId, s);

                    calcAllContributions();

                    sigZeroed = true;
                    break; // only remove the worst
                }
            }

            if(sigZeroed)
            {
                mRecentAllocPercents.clear(); // since allocations are likely to start moving again
                logStats();
            }

            // otherwise try again
            foundImprovements = findAdjustments();
            ++mIterations;

            if (!mIsValid)
                return false;

            if(mStagnantAllocChange)
                foundImprovements = false;

            if(targetSigZeroed)
                break;

            if (mIsFullyAllocated)
            {
                clearLowContribSigs();
                break;
            }

            if(mIterations >= MAX_TEST_ITERATIONS * 0.5)
            {
                if(!mLogVerboseOverride)
                {
                    // turn on logging to see what's happening
                    LOGGER.warn("sample({}) close to max iterations({}) reached", mSampleId, mIterations);
                    mLogVerboseOverride = true;
                }

                if(mIterations >= MAX_TEST_ITERATIONS)
                {
                    LOGGER.warn("sample({}) max iterations({}) reached", mSampleId, mIterations);
                    logStats();
                    break;
                }
            }
        }

        updateStats();
        return mIsValid;
    }

    private boolean findAdjustments()
    {
        List<Integer> exhaustedBuckets = getExhaustedBuckets(0.9);

        if (exhaustedBuckets.isEmpty())
        {
            calcAllContributions();
            exhaustedBuckets = getExhaustedBuckets(0.9);

            if (exhaustedBuckets.isEmpty())
                return false;
        }

        double[] maxOtherSigContribGains = new double[mSigCount];
        double maxReducedSigLoss = 0;
        double maxNetGain = 0;
        int maxReducedSig = -1;

        // find the sig with the max net gain across all other sigs from being reduced
        for (int sig = 0; sig < mSigCount; ++sig)
        {
            if(mContribs[sig] == 0 || (sig == mRequiredSig && mContribs[sig] <= mRequiredSigMinContrib))
                continue;

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
            // work out how many multiples of this combination of reductions and increases can be made before
            // either the reduced sig is exhausted or another sigs hits a limit elsewhere
            int applyMultiple = calcAdjustmentMultiple(maxReducedSig, maxReducedSigLoss, maxOtherSigContribGains);

            if(applyMultiple > 1)
            {
                maxReducedSigLoss *= applyMultiple;
                vectorMultiply(maxOtherSigContribGains, applyMultiple);
            }

            applySigAdjustments(maxReducedSig, maxReducedSigLoss, maxOtherSigContribGains);

            if (!mIsValid)
                return false;

            logStats();

            if (mIsFullyAllocated && mLogVerbose)
                LOGGER.debug(String.format("sample(%d) allocPerc(%.3f) exceeds target", mSampleId, mCurrentAllocPerc));

            return true;
        }

        return false;
    }

    private double testSigReduction(int rs, List<Integer> exhaustedBuckets, double[] otherSigContribGains)
    {
        double[][] sigData = mSigs.getData();

        double sig1ContribLoss = 0;
        double[] maxOtherSigContribs = new double[mSigCount];

        copyVector(mCurrentCounts, mReducedRefSigCounts);

        boolean initialTest = exhaustedBuckets.isEmpty();

        if(!initialTest)
        {
            // search through the exhausted buckets to find the one which delivers the largest gain in the reducing sig
            double minSig1ContribLoss = 0;
            for (Integer eb : exhaustedBuckets)
            {
                if (sigData[eb][rs] == 0)
                    continue;

                // calc the minimum the contribution loss in the exhausted bucket for the sig being reduced
                double contribLoss = mMinContribChange / sigData[eb][rs];

                if (minSig1ContribLoss == 0 || contribLoss < minSig1ContribLoss)
                {
                    minSig1ContribLoss = contribLoss;
                }
            }

            if (minSig1ContribLoss == 0)
                return 0;

            sig1ContribLoss = min(minSig1ContribLoss, mContribs[rs]); // cannot reduce past zero

            if(rs == mRequiredSig && mContribs[rs] - sig1ContribLoss < mRequiredSigMinContrib)
            {
                sig1ContribLoss = mContribs[rs] - mRequiredSigMinContrib;
            }

            // remove this sig across the board from a 1-lot contrib to this exhausted bucket
            for (int b = 0; b < mBucketCount; ++b)
            {
                mReducedRefSigCounts[b] -= sig1ContribLoss * sigData[b][rs];

                if (lessThan(mReducedRefSigCounts[b], 0))
                {
                    mIsValid = false;
                    LOGGER.error(String.format("bucket(%d) currentCount(%.1f) reduced below zero from sig(%d contrib=%.1f)",
                            b, mReducedRefSigCounts[b], mSigIds[rs], sig1ContribLoss));
                    return 0;
                }
            }
        }

        // now look at the potential gain to all the sigs in each bucket, having had the reduction sig's contribution removed
        for (int s2 = 0; s2 < mSigCount; ++s2)
        {
            if (rs == s2)
                continue;

            double minAlloc = calcSigContribution(s2, mReducedRefSigCounts);

            if (minAlloc > 0)
            {
                maxOtherSigContribs[s2] = minAlloc;
            }
        }

        double totalOtherSigsGain = sumVector(maxOtherSigContribs);

        if (totalOtherSigsGain < sig1ContribLoss * 0.5)
            return 0;

        // test out the proposed change to find the max that can be applied
        double[] testOtherSigContribs = new double[mSigCount];
        double maxOtherSigsGain = 0;

        // runs 0 and 1, sort ascending then descending with total allocation top-down each time
        // runs 2 and 3, sort ascending then descending with partial allocation top-down each time
        // runs 4 and 5, sort descending with BG first, with full then partial allocation
        for (int i = 0; i < 6; ++i)
        {
            copyVector(mReducedRefSigCounts, mTestRefSigCounts);

            // sort descending except on first iteration
            List<Integer> otherSigContribIndices = getSortedVectorIndices(maxOtherSigContribs, (i == 0 || i == 2));

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

                    double newAlloc = calcSigContribution(s2, mTestRefSigCounts);

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

                        if (greaterThan(mTestRefSigCounts[b] + newCount, mCounts[b]))
                        {
                            testOtherSigContribs[s2] = 0; // suggests an error in the calcSigContrib function
                            break;
                        }

                        mTestRefSigCounts[b] += newCount;
                    }
                }

                ++iterations;
            }

            // finally factor in potentially being able to add back in some portion of the sig that was reduced
            if(!initialTest)
            {
                double newAlloc = calcSigContribution(rs, mTestRefSigCounts);

                if (newAlloc > 0)
                {
                    testOtherSigContribs[rs] += newAlloc;
                }
            }

            double sigContribsTotal = sumVector(testOtherSigContribs);

            // don't take potential contributions which exclude the required sig
            if(initialTest && mRequiredSig >= 0 && testOtherSigContribs[mRequiredSig] == 0)
                continue;

            if (sigContribsTotal > maxOtherSigsGain)
            {
                // take the top allocation combination
                maxOtherSigsGain = sigContribsTotal;
                copyVector(testOtherSigContribs, otherSigContribGains);
            }
        }

        if (maxOtherSigsGain < sig1ContribLoss)
            return 0;

        // check for minisule gains and loss
        if(maxOtherSigsGain < 0.001)
            return 0;

        return sig1ContribLoss;
    }

    private int calcAdjustmentMultiple(int rs, double sigContribLoss, double[] otherSigContribGains)
    {
        // given the proposed reduction of a sig's contribution (ie sigContribLoss), for expediency determine if this
        // reduction can be made multiple times in one go and if so to what extent (the return integer)
        copyVector(mCurrentCounts, mTestCounts);
        double[] testContribs = new double[mSigCount];
        copyVector(mContribs, testContribs);

        int applyMultiple = 0;
        int maxMultiples = 100;
        boolean appliedOk = true;

        final double[][] sigData = mSigs.getData();

        while(appliedOk)
        {
            // reduce the main sig
            for (int b = 0; b < mBucketCount; ++b)
            {
                double newCount = sigContribLoss * sigData[b][rs];

                if (lessThan(mTestCounts[b] - newCount, 0))
                {
                    appliedOk = false;
                    break;
                }

                mTestCounts[b] -= newCount;
            }

            if(!appliedOk)
                break;

            if(lessThan(testContribs[rs] - sigContribLoss, 0))
                break;

            if(rs == mRequiredSig && testContribs[rs] - sigContribLoss < mRequiredSigMinContrib)
                break;

            testContribs[rs] -= sigContribLoss;

            // increase all the others
            for(int sig = 0; sig < mSigCount; ++sig)
            {
                if(otherSigContribGains[sig] == 0)
                    continue;

                for (int b = 0; b < mBucketCount; ++b)
                {
                    double newCount = otherSigContribGains[sig] * sigData[b][sig];

                    if (greaterThan(mTestCounts[b] + newCount, mCounts[b]))
                    {
                        appliedOk = false;
                        break;
                    }

                    mTestCounts[b] += newCount;
                }

                testContribs[sig] += otherSigContribGains[sig];
            }

            if(!appliedOk)
                break;

            if(greaterThan(sumVector(testContribs), mCountsTotal))
                break;

            ++applyMultiple;

            if(applyMultiple >= maxMultiples)
                break;
        }

        return applyMultiple;
    }

    private void applySigAdjustments(int rs, double sigContribLoss, double[] otherSigContribGains)
    {
        if(rs >= 0)
        {
            if (log())
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
            if (otherSigContribGains[s2] == 0)
                break;

            if (log())
            {
                LOGGER.debug(String.format("apply sig(%d) contrib(%s -> %s gain=%s)",
                        mSigIds[s2], doubleToStr(mContribs[s2]), doubleToStr(mContribs[s2] + otherSigContribGains[s2]), doubleToStr(otherSigContribGains[s2])));
            }

            applyContribution(s2, otherSigContribGains[s2]);
        }
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
        double minAlloc = 0;
        double[][] sigData = mSigs.getData();
        double contribLimit = mCountsTotal - mContribTotal;

        for (int b = 0; b < mBucketCount; ++b)
        {
            if (sigData[b][sig] == 0)
                continue;

            if (mCounts[b] > 0 && currentCounts[b] >= mCounts[b])
            {
                minAlloc = 0;
                break;
            }

            double sigRatio = sigData[b][sig];
            double alloc = (mCounts[b] - currentCounts[b]) / sigRatio;

            alloc = min(alloc, mCountsTotal);
            // alloc = min(alloc, contribLimit);

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
        if (greaterThan(mContribTotal + newContrib, mCountsTotal))
        {
            if(log())
            {
                LOGGER.debug(String.format("sample(%d) lowering invalid contrib addition(%.1f + %.1f) vs countsLimit(%.1f)",
                        mSampleId, mContribTotal, newContrib, mCountsTotal));
            }

            // no longer invalid since the capping of noise can cause this
            // mIsValid = false;
            // return;

            newContrib = mCountsTotal - mContribTotal;

            if(newContrib == 0)
                return;
        }

        if (lessThan(mContribTotal + newContrib, 0))
        {
            LOGGER.error(String.format("sample(%d) invalid contrib addition(%.1f + %.1f)",
                    mSampleId, mContribTotal, newContrib, mCountsTotal));
            mIsValid = false;
            return;
        }

        if (lessThan(mContribs[sig] + newContrib, 0))
        {
            LOGGER.error(String.format("sample(%d) invalid sig(%d) contrib reduction(%.1f + %.1f)",
                    mSampleId, mSigIds[sig], mContribs[sig], newContrib));
            mIsValid = false;
            return;
        }

        mContribs[sig] += newContrib;
        mContribTotal += newContrib;

        final double[][] sigData = mSigs.getData();

        for (int b = 0; b < mBucketCount; ++b)
        {
            if(sigData[b][sig] == 0)
                continue;

            double newCount = newContrib * sigData[b][sig];

            if (greaterThan(mCurrentCounts[b] + newCount, mCounts[b]))
            {
                LOGGER.error(String.format("sample(%d) sig(%d contrib=%f) count(%f) + newCount(%f) exceeds maxCount(%f) diff(%f)",
                        mSampleId, mSigIds[sig], mContribs[sig], mCurrentCounts[b], newCount, mCounts[b],
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
            if(mContribs[sig] > 1)
            {
                LOGGER.error("sample({}) sig({}) previously zeroed", mSampleId, mSigIds[sig]);
                mIsValid = false;
            }
            return;
        }

        mZeroedSigs.add(sig);

        double sigContrib = mContribs[sig];

        double[][] sigData = mSigs.getData();

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

            sigData[b][sig] = 0;
        }

        mContribTotal -= mContribs[sig];

        mContribs[sig] = 0;

        if (log())
            LOGGER.debug(String.format("sample(%d) sig(%d) contrib(%s) zeroed", mSampleId, mSigIds[sig], doubleToStr(sigContrib)));
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

        /* no longer applicable since mCountsTotal is capped by max noise
        if (abs(mCountsTotal - mContribTotal - mResiduals) >= 1)
        {
            LOGGER.error(String.format("sample(%d) totalCount(%.1f) less contribTotal(%.1f) != residuals(%.1f))",
                    mSampleId, mCountsTotal, mContribTotal, mResiduals));
            mIsValid = false;
        }
        */

        List<Integer> sortedContribIndices = getSortedVectorIndices(mContribs, false);

        mHasLowContribSigs = false;

        initVector(mCurrentCountsNoNoise, 0);

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
                    mCurrentCountsNoNoise[i] = min(mCurrentCountsNoNoise[i] + (mContribs[s] * sigData[i][s]), mRawCounts[i]);
                }
            }

            if (s == mRequiredSig)
                continue;

            if (!aboveMinReqContrib(mContribs[s]))
                mHasLowContribSigs = true;
        }

        mCurrentAllocTotal = sumVector(mCurrentCountsNoNoise);
        mCurrentAllocPerc = min(mCurrentAllocTotal / mRawCountsTotal, 1);

        // if say the target is 99%, then residuals must also be sufficiently small
        mIsFullyAllocated = mCurrentAllocPerc >= mTargetAllocPercent && mResiduals/mRawCountsTotal <= mTargetResidualsPercent;

        if(mInitAllocPerc == 0)
        {
            mInitAllocPerc = mCurrentAllocPerc;
            copyVector(mContribs, mInitContribs);
        }

        mRecentAllocPercents.add(mCurrentAllocPerc);

        if(mRecentAllocPercents.size() >= MAX_NO_IMPROVE_COUNT)
        {
            double maxVal = 0;
            double minVal = 1;
            for(Double allocPerc : mRecentAllocPercents)
            {
                maxVal = max(allocPerc, maxVal);
                minVal = min(allocPerc, minVal);
            }

            double range = maxVal - minVal;
            double recentMove = mCurrentAllocPerc - mRecentAllocPercents.get(0);

            if(range < 0.01 || recentMove < 0.01)
            {
                mStagnantAllocChange = true;

                if(mTargetSig == -1)
                {
                    LOGGER.debug(String.format("sample(%d) stagnant recent %d moves: range(%.1f min=%.1f max=%.1f start=%.1f end=%.1f)",
                            mSampleId, mRecentAllocPercents.size(), range, minVal, maxVal, mRecentAllocPercents.get(0), mCurrentAllocPerc));

                    mLogVerboseOverride = true;
                }
            }

            mRecentAllocPercents.remove(0);
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

        mMinContribChange = max(mRawCountsTotal * requiredChgPerc, absMinChange);
        */

        mMinContribChange = max(mRawCountsTotal * MIN_COUNT_CHG_PERC, MIN_COUNT_CHG);
    }

    private boolean aboveMinReqContrib(double contrib)
    {
        return (contrib / mRawCountsTotal >= mMinContribPercent) && contrib >= mMinContribCount;
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

    private boolean log() { return mLogVerbose || mLogVerboseOverride; }

    private void logStats()
    {
        recalcStats();

        if (!log())
            return;

        double underAllocation = max(mRawCountsTotal - mContribTotal, 0);

        LOGGER.debug(String.format("sample(%d) totalCount(%s maxNs=%s) allocated(%s wNs=%.3f noNs=%.3f init=%.3f) underAlloc(%s perc=%.3f) res(%s perc=%.3f) data(mc=%.1f ls=%s it=%d)",
                mSampleId, doubleToStr(mRawCountsTotal), doubleToStr(mMaxNoiseTotal), doubleToStr(mContribTotal), mContribTotal / mRawCountsTotal, mCurrentAllocPerc, mInitAllocPerc,
                doubleToStr(underAllocation), underAllocation/mRawCountsTotal, doubleToStr(mResiduals), mResiduals/mRawCountsTotal,
                mMinContribChange, mHasLowContribSigs, mIterations));

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

            contribStr += String.format("%d = %s perc=%.3f", mSigIds[s], doubleToStr(mContribs[s]), min(mContribs[s] / mRawCountsTotal, 1));
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

    private void updateStats()
    {
        if(!mIsValid)
            return;

        mAvgIterations = (mInstances * mAvgIterations + mIterations) / (double)(mInstances + 1);
        double percImprove = max(mCurrentAllocPerc - mInitAllocPerc, 0);
        mAvgPercImprove = (mInstances * mAvgPercImprove + percImprove) / (mInstances + 1);
        ++mInstances;
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
            mResiduals += abs(mRawCounts[i] - mCurrentCounts[i]);
        }
    }

}