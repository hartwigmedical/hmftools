package com.hartwig.hmftools.sigs.buckets;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.sigs.DataUtils.doubleToStr;
import static com.hartwig.hmftools.common.sigs.DataUtils.doublesEqual;
import static com.hartwig.hmftools.common.sigs.DataUtils.greaterThan;
import static com.hartwig.hmftools.common.sigs.DataUtils.lessThan;
import static com.hartwig.hmftools.common.sigs.DataUtils.sizeToStr;
import static com.hartwig.hmftools.common.utils.VectorUtils.copyVector;
import static com.hartwig.hmftools.common.utils.VectorUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.common.utils.VectorUtils.initVector;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVectors;
import static com.hartwig.hmftools.common.utils.VectorUtils.vectorMultiply;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SampleSigContribOptimiser
{
    private SampleData mSample;
    private final List<List<Integer>> mBucketIdsCollection; // the signature bucket ratios
    private final List<double[]> mSigAllocCounts;
    private final List<double[]> mTestSigNewCounts;
    private final List<double[]> mOtherSigNewCounts;
    private final List<double[]> mMaxOtherSigNewCounts;
    private final double[] mMaxRefSigReductionCounts; // reduction to counts in the sig being tested
    private final double[] mCurrentAllocCounts;
    private final double[] mCurrentAllocNoise;

    private int mSampleId;
    private boolean mIsValid;

    private final double[] mRawCounts; // actual sample counts
    private final double[] mCountsNoise; // noise around the counts
    private final List<double[]> mRatiosCollection; // the signature bucket ratios
    private int[] mSigIds;

    private double[] mInitContribs;
    private double mMinContribPercent; // each sig's final contrib must be above this level as percent of total variants
    private int mMinContribCount; // each sig's final contrib must be above this level in absolute terms
    private double mTargetAllocPercent; // target total allocation to exit the fit
    private double mTargetResidualsPercent;
    private int mTargetSig; // to check if a specific sig remains above the require contribution percent
    private int mRequiredSig; // keep this specific sig even if it falls below the require contribution percent
    private double mRequiredSigMinContrib;
    private final List<Integer> mZeroedSigs; // the signature bucket ratios

    private final double[] mCounts; // sample counts, optionally including noise
    private double[] mContribs;
    private double mContribTotal;

    private final int mBucketCount;
    private int mSigCount;
    private double mRawCountsTotal; // total of the actual variant counts
    private double mResiduals; // standard meaning - absolute value of diff between actual and allocated (including excess)
    private double mCurrentAllocTotal; // sum of contributions above the require percent and capped at actual counts, not noise
    private double mCurrentAllocPerc;
    private double mInitAllocPerc;
    private boolean mHasLowContribSigs;
    private boolean mIsFullyAllocated;

    // calc state to avoid repeated memory alloc
    private double[] mReducedRefSigCounts; // reduction to counts in the sig being tested
    private final double[] mApplyMultipleCounts; // counts to test applying a multiple of the reduction/increase routine

    // diagnostics and stats
    private int mIterations;
    private int mInstances;
    private double mAvgIterations;
    private double mAvgPercImprove;
    private final List<Double> mRecentAllocPercents;
    private boolean mStagnantAllocChange;

    private double mMinContribChange;
    private double mMinContribChangePercent;
    private double mMinContribChangeValue;

    private boolean mLogVerbose;
    private boolean mLogVerboseOverride;

    // constants and control config
    private static double MIN_COUNT_CHG_PERC = 0.001;
    private static double MIN_COUNT_CHG = 1;
    private static int MAX_TEST_ITERATIONS = 100;
    private static int MAX_NO_IMPROVE_COUNT = 10;
    private static double REQUIRED_SIG_MIN_PERCENT = 0.001;

    private static final Logger LOGGER = LogManager.getLogger(SampleSigContribOptimiser.class);

    public SampleSigContribOptimiser(int bucketCount, boolean logVerbose, double targetAllocPercent)
    {
        mSample = null;
        mBucketCount = bucketCount;

        mLogVerbose = logVerbose;
        mMinContribPercent = 0.001;
        mTargetAllocPercent = targetAllocPercent;
        mTargetResidualsPercent = (1 - targetAllocPercent) * 5;

        mMinContribChangePercent = MIN_COUNT_CHG_PERC;
        mMinContribChangeValue = MIN_COUNT_CHG;

        mTargetSig = -1;
        mRequiredSig = -1;
        mRequiredSigMinContrib = 0;
        mZeroedSigs = Lists.newArrayList();
        mRatiosCollection = Lists.newArrayList();
        mRecentAllocPercents = Lists.newArrayList();

        mRawCounts = new double[mBucketCount];
        mCounts = new double[mBucketCount];
        mCountsNoise = new double[mBucketCount];

        mReducedRefSigCounts = new double[mBucketCount];
        mApplyMultipleCounts = new double[mBucketCount];

        mCurrentAllocCounts = new double[mBucketCount];
        mCurrentAllocNoise = new double[mBucketCount];
        mMaxRefSigReductionCounts = new double[mBucketCount];

        mSigAllocCounts = Lists.newArrayList();
        mOtherSigNewCounts = Lists.newArrayList();
        mTestSigNewCounts = Lists.newArrayList();
        mMaxOtherSigNewCounts = Lists.newArrayList();
        mBucketIdsCollection = Lists.newArrayList();

        mInstances = 0;
        mAvgIterations = 0;
        mAvgPercImprove = 0;
    }

    public void initialise(final SampleData sample, final List<double[]> ratiosCollection, double minSigPercent, int minAllocCount)
    {
        mSample = sample;
        mSampleId = sample.Id;
        mMinContribPercent = minSigPercent;
        mMinContribCount = minAllocCount;
        mMinContribChange = 0;

        mSigAllocCounts.clear();
        mOtherSigNewCounts.clear();
        mTestSigNewCounts.clear();
        mMaxOtherSigNewCounts.clear();
        mBucketIdsCollection.clear();

        mSigCount = ratiosCollection.size();
        mSigIds = new int[mSigCount];
        for(int i = 0; i < mSigCount; ++i)
        {
            mSigIds[i] = i;
            mSigAllocCounts.add(new double[mBucketCount]);
            mOtherSigNewCounts.add(new double[mBucketCount]);
            mTestSigNewCounts.add(new double[mBucketCount]);
            mMaxOtherSigNewCounts.add(new double[mBucketCount]);
        }

        initVector(mCurrentAllocCounts, 0);
        initVector(mCurrentAllocNoise, 0);
        initVector(mMaxRefSigReductionCounts, 0);
        initVector(mReducedRefSigCounts, 0);
        initVector(mApplyMultipleCounts, 0);
        initVector(mCounts, 0);

        if(sample.usingElevatedForAllocation())
            copyVector(sample.getElevatedBucketCounts(), mRawCounts);
        else
            copyVector(sample.getBucketCounts(), mRawCounts);

        copyVector(mRawCounts, mCounts);
        copyVector(sample.getNoiseCounts(), mCountsNoise);

        sumVectors(mCountsNoise, mCounts);

        mTargetSig = -1;
        mRequiredSig = -1;
        mRequiredSigMinContrib = 0;
        mZeroedSigs.clear();

        mRatiosCollection.clear();
        mRatiosCollection.addAll(ratiosCollection);

        mContribs = new double[mSigCount];
        mInitContribs = new double[mSigCount];
        mContribTotal = 0;

        mRawCountsTotal = sumVector(mRawCounts);

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

        if (mContribs.length != ratiosCollection.size())
        {
            mIsValid = false;
            return;
        }

        for (int sig = 0; sig < mSigCount; ++sig)
        {
            double[] sigRatios = ratiosCollection.get(sig);

            List<Integer> bucketIds = Lists.newArrayList();
            mBucketIdsCollection.add(bucketIds);

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
                if(sigRatios[b] > 0)
                    bucketIds.add(b);
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

    public void setMinContribChange(double value, double percent)
    {
        mMinContribChangePercent = percent;
        mMinContribChangeValue = value;
    }

    public final double[] getContribs() { return mContribs; }
    public List<double[]> getSigAllocCounts() { return mSigAllocCounts; }
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

            if(mInitAllocPerc == 0)
                logStats();

            exhaustedBuckets = getExhaustedBuckets(0.9);

            if (exhaustedBuckets.isEmpty())
                return false;
        }

        double[] maxOtherSigContribGains = new double[mSigCount];
        double[] otherSigContribGains = new double[mSigCount];
        double maxReducedSigLoss = 0;
        double maxNetGain = 0;
        int maxReducedSig = -1;

        copyVector(mSample.getAllocBucketCounts(), mCurrentAllocCounts);
        copyVector(mSample.getAllocNoiseCounts(), mCurrentAllocNoise);
        initVector(mMaxRefSigReductionCounts, 0);

        for(int i = 0; i < mSigCount; ++i)
        {
            double[] maxOtherCounts = mMaxOtherSigNewCounts.get(i);
            initVector(maxOtherCounts, 0);
        }

        // find the sig with the max net gain across all other sigs from being reduced
        for (int sig = 0; sig < mSigCount; ++sig)
        {
            if(mContribs[sig] == 0 || (sig == mRequiredSig && mContribs[sig] <= mRequiredSigMinContrib))
                continue;

            initVector(otherSigContribGains, 0);
            double reducedSigContribLoss = testSigReduction(sig, exhaustedBuckets, otherSigContribGains);

            // restore the counts from this testing routine
            mSample.restoreCounts(mCurrentAllocCounts, mCurrentAllocNoise);

            double netGain = sumVector(otherSigContribGains) - reducedSigContribLoss;

            if (reducedSigContribLoss > 0 && netGain > maxNetGain)
            {
                maxNetGain = netGain;
                maxReducedSig = sig;
                maxReducedSigLoss = reducedSigContribLoss;
                copyVector(otherSigContribGains, maxOtherSigContribGains);
                copyVector(mReducedRefSigCounts, mMaxRefSigReductionCounts);

                for(int i = 0; i < mSigCount; ++i)
                {
                    double[] maxOtherCounts = mMaxOtherSigNewCounts.get(i);
                    double[] otherCounts = mOtherSigNewCounts.get(i);
                    copyVector(otherCounts, maxOtherCounts);
                }
            }
        }

        if (maxReducedSig >= 0)
        {
            // work out how many multiples of this combination of reductions and increases can be made before
            // either the reduced sig is exhausted or another sigs hits a limit elsewhere
            int applyMultiple = calcAdjustmentMultiple(maxReducedSig, maxReducedSigLoss, maxOtherSigContribGains);

            if(applyMultiple > 1)
            {
                vectorMultiply(mMaxRefSigReductionCounts, applyMultiple);

                maxReducedSigLoss *= applyMultiple;
                vectorMultiply(maxOtherSigContribGains, applyMultiple);

                for(int i = 0; i < mSigCount; ++i)
                {
                    double[] maxOtherCounts = mMaxOtherSigNewCounts.get(i);
                    vectorMultiply(maxOtherCounts, applyMultiple);
                }
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

    private static double REDUCED_ALLOC_PERCENT = 0.75;

    private double testSigReduction(int reductSig, List<Integer> exhaustedBuckets, double[] otherSigContribGains)
    {
        // reduce the reference sig (rs) by enough to expose unallocated counts for the othet sigs
        // return and set:
        // - reductions to the ref sig's counts
        // - the sum of the gains to the other sigs (than the ref)
        // - the set of new contributions
        // - the other sig's new counts
        double sig1ContribLoss = 0;
        double[] potentialOtherSigContribs = new double[mSigCount];

        initVector(mReducedRefSigCounts, 0);

        double[] allocCounts = new double[mBucketCount]; // will function as a spare array for various purposes

        boolean initialTest = exhaustedBuckets.isEmpty();

        if(!initialTest)
        {
            final double[] rsRatios = mRatiosCollection.get(reductSig);

            // search through the exhausted buckets to find the one which delivers the largest gain in the reducing sig
            double minSig1ContribLoss = 0;
            for (Integer eb : exhaustedBuckets)
            {
                if (rsRatios[eb] == 0)
                    continue;

                // calc the minimum the contribution loss in the exhausted bucket for the sig being reduced
                double contribLoss = mMinContribChange / rsRatios[eb];

                if (minSig1ContribLoss == 0 || contribLoss < minSig1ContribLoss)
                {
                    minSig1ContribLoss = contribLoss;
                }
            }

            if (minSig1ContribLoss == 0)
                return 0;

            sig1ContribLoss = min(minSig1ContribLoss, mContribs[reductSig]); // cannot reduce past zero

            if(reductSig == mRequiredSig && mContribs[reductSig] - sig1ContribLoss < mRequiredSigMinContrib)
            {
                sig1ContribLoss = mContribs[reductSig] - mRequiredSigMinContrib;
            }

            // remove this sig across the board from a 1-lot contrib to this exhausted bucket
            final double[] sigAllocCounts = mSigAllocCounts.get(reductSig);
            double actualSigLoss = 0;

            for (int b = 0; b < mBucketCount; ++b)
            {
                mReducedRefSigCounts[b] = max(-sig1ContribLoss * rsRatios[b], -sigAllocCounts[b]);
                actualSigLoss += mReducedRefSigCounts[b];
            }

            sig1ContribLoss = -actualSigLoss;
        }

        // now look at the potential gain to all the sigs in each bucket, having had the reduction sig's contribution removed
        mSample.reduceAllocCounts(mReducedRefSigCounts);

        for (int s2 = 0; s2 < mSigCount; ++s2)
        {
            if (reductSig == s2 || mZeroedSigs.contains(s2))
                continue;

            double potentialAlloc = mSample.getPotentialUnallocCounts(mRatiosCollection.get(s2), mBucketIdsCollection.get(s2),
                    null, allocCounts);

            if (potentialAlloc > 0)
            {
                potentialOtherSigContribs[s2] = potentialAlloc;
            }
        }

        // ref sig's counts will remain reduced for the next phase
        double totalOtherSigsGain = sumVector(potentialOtherSigContribs);

        if (totalOtherSigsGain < sig1ContribLoss * 0.5)
            return 0;

        // test out the proposed change to find the max that can be applied
        double[] testOtherSigContribs = new double[mSigCount];
        double maxOtherSigsGain = 0;
        int maxIterationIndex = 0;

        // runs 0 and 1, sort ascending then descending with total allocation top-down each time
        // runs 2 and 3, sort ascending then descending with partial allocation top-down each time
        // runs 4 and 5, sort descending with any required sig put first, with full then partial allocation
        for (int i = 0; i < 6; ++i)
        {
            if(i > 0)
            {
                // revert back to the start point, where only the ref sig's counts have been reduced
                mSample.restoreCounts(mCurrentAllocCounts, mCurrentAllocNoise);
                mSample.reduceAllocCounts(mReducedRefSigCounts);

                for (int j = 0; j < mSigCount; ++j)
                {
                    double[] testOtherCounts = mTestSigNewCounts.get(j);
                    initVector(testOtherCounts, 0);
                }
            }

            // sort descending except on first iteration
            List<Integer> otherSigContribIndices = getSortedVectorIndices(potentialOtherSigContribs, (i == 0 || i == 2));

            // try the likely order from the original discovery fit process
            if(i >= 4)
            {
                if(!initialTest || mRequiredSig == -1)
                    break;

                if(otherSigContribIndices.get(0) == mRequiredSig)
                    break;

                // put the required sig first
                otherSigContribIndices.add(0, mRequiredSig);
            }

            if (i > 0)
                initVector(testOtherSigContribs, 0);

            boolean usingReducedAllocation = (i == 2 || i == 3 || i == 5);
            double allocPerc = usingReducedAllocation ? REDUCED_ALLOC_PERCENT : 1.0;
            boolean foundAdjusts = true;
            int iterations = 0;
            int maxIteration = 5;

            // reduced allocation is leading to inconsistent results, since in this routine sigs are allocated repeatedly, rather than once each
            // using their cumulative contribution
            if(usingReducedAllocation)
                continue;

            while (foundAdjusts && iterations < maxIteration)
            {
                foundAdjusts = false;

                for (Integer s2 : otherSigContribIndices)
                {
                    if (s2 == reductSig || potentialOtherSigContribs[s2] == 0)
                        continue;

                    double potentialAlloc = mSample.getPotentialUnallocCounts(mRatiosCollection.get(s2), mBucketIdsCollection.get(s2),
                            null, allocCounts);

                    if (potentialAlloc <= 0)
                        continue;

                    if (iterations < maxIteration - 1 && allocPerc < 1)
                    {
                        potentialAlloc *= allocPerc;
                        vectorMultiply(allocCounts, allocPerc);
                    }

                    // and then apply them
                    double actualAlloc = mSample.allocateBucketCounts(allocCounts);

                    if(!doublesEqual(actualAlloc, potentialAlloc, 0.1))
                    {
                        LOGGER.warn(String.format("sample(%d) potentialAlloc(%.1f) != actualAlloc(%.1f)", mSampleId, potentialAlloc, actualAlloc));
                        return 0;
                    }

                    double[] otherSigNewCounts = mTestSigNewCounts.get(s2);
                    testOtherSigContribs[s2] += actualAlloc; // take adjustment if required
                    sumVectors(allocCounts, otherSigNewCounts);

                    foundAdjusts = true;
                }

                if(!usingReducedAllocation)
                    break;

                ++iterations;
            }

            // finally factor in potentially being able to add back in some portion of the sig that was reduced
            if(!initialTest)
            {
                double potentialAlloc = mSample.getPotentialUnallocCounts(mRatiosCollection.get(reductSig), mBucketIdsCollection.get(reductSig),
                        null, allocCounts);

                if (potentialAlloc > 0)
                {
                    double actualAlloc = mSample.allocateBucketCounts(allocCounts);

                    if(!doublesEqual(actualAlloc, potentialAlloc, 0.1))
                    {
                        LOGGER.warn(String.format("sample(%d) potentialAlloc(%.1f) != actualAlloc(%.1f)", mSampleId, potentialAlloc, actualAlloc));
                        return 0;
                    }

                    testOtherSigContribs[reductSig] += actualAlloc;
                    double[] otherSigNewCounts = mTestSigNewCounts.get(reductSig);
                    sumVectors(allocCounts, otherSigNewCounts);
                }
            }

            double sigContribsTotal = sumVector(testOtherSigContribs);

            // don't take potential contributions which exclude the required sig
            if(initialTest && mRequiredSig >= 0 && testOtherSigContribs[mRequiredSig] == 0)
                continue;

            if (sigContribsTotal > maxOtherSigsGain)
            {
                // take the top allocation combination
                maxIterationIndex = i;
                maxOtherSigsGain = sigContribsTotal;
                copyVector(testOtherSigContribs, otherSigContribGains);

                for(int j = 0; j < mSigCount; ++j)
                {
                    double[] testOtherCounts = mTestSigNewCounts.get(j);
                    double[] otherCounts = mOtherSigNewCounts.get(j);
                    copyVector(testOtherCounts, otherCounts);
                }
            }
        }

        if (maxOtherSigsGain < sig1ContribLoss)
            return 0;

        // check for minisule gains and loss
        if(maxOtherSigsGain < 0.001)
            return 0;

        if(log())
        {
            LOGGER.debug("sample({}) max gain({}) for refSigLoss({}) at iterationType({})",
                    mSampleId, sizeToStr(maxOtherSigsGain), sizeToStr(sig1ContribLoss), maxIterationIndex);
        }

        return sig1ContribLoss;
    }

    private int calcAdjustmentMultiple(int reductSig, double sigContribLoss, double[] otherSigContribGains)
    {
        // given the proposed reduction of a sig's contribution (ie sigContribLoss), for expediency determine if this
        // reduction can be made multiple times in one go and if so to what extent (the return integer)
        copyVector(mCurrentAllocCounts, mApplyMultipleCounts);
        double[] testContribs = new double[mSigCount];
        copyVector(mContribs, testContribs);

        int applyMultiple = 0;
        int maxMultiples = 100;
        boolean appliedOk = true;

        while(appliedOk)
        {
            // reduce the main sig
            final double[] rsRatios = mRatiosCollection.get(reductSig);

            for (int b = 0; b < mBucketCount; ++b)
            {
                double newCount = sigContribLoss * rsRatios[b];

                if (lessThan(mApplyMultipleCounts[b] - newCount, 0))
                {
                    appliedOk = false;
                    break;
                }

                mApplyMultipleCounts[b] -= newCount;
            }

            if(!appliedOk)
                break;

            if(lessThan(testContribs[reductSig] - sigContribLoss, 0))
                break;

            if(reductSig == mRequiredSig && testContribs[reductSig] - sigContribLoss < mRequiredSigMinContrib)
                break;

            testContribs[reductSig] -= sigContribLoss;

            // increase all the others
            for(int sig = 0; sig < mSigCount; ++sig)
            {
                if(otherSigContribGains[sig] == 0)
                    continue;

                final double[] sigRatios = mRatiosCollection.get(sig);

                for (int b = 0; b < mBucketCount; ++b)
                {
                    double newCount = otherSigContribGains[sig] * sigRatios[b];

                    if (greaterThan(mApplyMultipleCounts[b] + newCount, mCounts[b]))
                    {
                        appliedOk = false;
                        break;
                    }

                    mApplyMultipleCounts[b] += newCount;
                }

                testContribs[sig] += otherSigContribGains[sig];
            }

            if(!appliedOk)
                break;

            ++applyMultiple;

            if(applyMultiple >= maxMultiples)
                break;
        }

        return applyMultiple;
    }

    private void applySigAdjustments(int reductSig, double refSigContribLoss, double[] otherSigContribGains)
    {
        if(reductSig >= 0)
        {
            if (log())
            {
                double totalActualOtherGain = sumVector(otherSigContribGains);

                LOGGER.debug(String.format("reduceSig(%s cur=%s loss=%s) for other sigs gain(%s)",
                        mSigIds[reductSig], doubleToStr(mContribs[reductSig]), doubleToStr(refSigContribLoss), doubleToStr(totalActualOtherGain)));
            }

            applyContribution(reductSig, mMaxRefSigReductionCounts, -refSigContribLoss);
        }

        final List<Integer> otherSigContribIndices = getSortedVectorIndices(otherSigContribGains, false);

        double[] otherSigCounts = new double[mBucketCount];

        for (Integer otherSig : otherSigContribIndices)
        {
            double otherSigContrib = otherSigContribGains[otherSig];
            if (otherSigContrib == 0)
                break;

            if (log())
            {
                LOGGER.debug(String.format("apply sig(%d) contrib(%s -> %s gain=%s)",
                        mSigIds[otherSig], doubleToStr(mContribs[otherSig]), doubleToStr(mContribs[otherSig] + otherSigContrib),
                        doubleToStr(otherSigContribGains[otherSig])));
            }

            final double[] sigDefn = mRatiosCollection.get(otherSig);

            for(int b = 0; b < mBucketCount; ++b)
            {
                otherSigCounts[b] = otherSigContrib * sigDefn[b];
            }

            applyContribution(otherSig, mMaxOtherSigNewCounts.get(otherSig), otherSigContribGains[otherSig]);
            // applyContribution(otherSig, otherSigCounts, otherSigContrib);
        }
    }

    private void calcAllContributions()
    {
        double[] allocCounts = new double[mBucketCount];
        for (int s = 0; s < mSigCount; ++s)
        {
            if (mZeroedSigs.contains(s))
                continue;

            double allocTotal = calcSigContribution(s, allocCounts);
            applyContribution(s, allocCounts, allocTotal);
        }
    }

    private double calcSigContribution(int sig, double[] allocCounts)
    {
        return mSample.getPotentialUnallocCounts(mRatiosCollection.get(sig), mBucketIdsCollection.get(sig), null, allocCounts);
    }

    private void applyContribution(int sig, double[] newCounts, double newCountsTotal)
    {
        if(newCountsTotal == 0)
            return;

        if(!doublesEqual(newCountsTotal, sumVector(newCounts), 0.1))
        {
            LOGGER.warn(String.format("sample(%d) about to apply inconsistent newCountsTotal(%.1f) != allocChange(%.1f)",
                    mSampleId, newCountsTotal, abs(sumVector(newCounts))));
            mIsValid = false;
        }

        double allocChange = 0;

        if(newCountsTotal > 0)
        {
            allocChange = mSample.allocateBucketCounts(newCounts);
        }
        else
        {
            allocChange = mSample.reduceAllocCounts(newCounts);
        }

        if(!doublesEqual(newCountsTotal, allocChange, 0.1))
        {
            LOGGER.warn(String.format("sample(%d) newCountsTotal(%.1f) != allocChange(%.1f)", mSampleId, newCountsTotal, allocChange));
            mIsValid = false;
        }

        double[] sigAllocCounts = mSigAllocCounts.get(sig);
        sumVectors(newCounts, sigAllocCounts);

        if(log())
        {
            for(int i = 0; i < sigAllocCounts.length; ++i)
            {
                if(sigAllocCounts[i] < 0)
                {
                    LOGGER.warn(String.format("sample(%d) bucket(%d) has negative count(%.1f)", mSampleId, i, sigAllocCounts[i]));
                }
            }
        }

        mContribs[sig] += allocChange;
        mContribTotal += allocChange;

        // update local state to match the sample
        copyVector(mSample.getAllocBucketCounts(), mCurrentAllocCounts);
        copyVector(mSample.getAllocNoiseCounts(), mCurrentAllocNoise);
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

        double[] sigCounts = new double[mBucketCount];
        copyVector(mSigAllocCounts.get(sig), sigCounts);
        vectorMultiply(sigCounts, -1);
        double reductionTotal = sumVector(sigCounts);

        applyContribution(sig, sigCounts, reductionTotal);

        mContribs[sig] = 0;

        if (log())
            LOGGER.debug(String.format("sample(%d) sig(%d) contrib(%s) zeroed", mSampleId, mSigIds[sig], doubleToStr(sigContrib)));
    }

    private void recalcStats()
    {
        // validate current state
        calcResiduals();

        mContribTotal = sumVector(mContribs);

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

        mCurrentAllocTotal = mSample.getAllocatedCount();
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
        // make the min contrib change a function of the current and target alloc
        // to have it move more aggressively when further out and then more fine-grained close to the target
        double targetRemaining = max(mTargetAllocPercent - mCurrentAllocPerc, 0) / mTargetAllocPercent;

        double upperPercent = 0.05;
        double lowerPercent = 0.001;
        double absMinChange = 0.3;
        double requiredChgPerc = lowerPercent + (upperPercent - lowerPercent) * targetRemaining;

        mMinContribChange = max(mRawCountsTotal * requiredChgPerc, absMinChange);

        // mMinContribChange = max(mRawCountsTotal * mMinContribChangePercent, mMinContribChangeValue);
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

        LOGGER.debug(String.format("sample(%d) totalCount(%s) allocated(%s wNs=%.3f noNs=%.3f init=%.3f) underAlloc(%s perc=%.3f) res(%s perc=%.3f) data(mc=%.1f ls=%s it=%d)",
                mSampleId, doubleToStr(mRawCountsTotal), doubleToStr(mContribTotal), mContribTotal / mRawCountsTotal, mCurrentAllocPerc, mInitAllocPerc,
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
            if (mRawCounts[b] > 0 && mSample.getAllocBucketCounts()[b] >= percFull * mRawCounts[b])
                exhaustedBuckets.add(b);
        }

        return exhaustedBuckets;
    }

    private void calcResiduals()
    {
        mResiduals = 0;

        for (int i = 0; i < mBucketCount; ++i)
        {
            mResiduals += abs(mRawCounts[i] - (mCurrentAllocCounts[i] + mCurrentAllocNoise[i]));
        }
    }

}