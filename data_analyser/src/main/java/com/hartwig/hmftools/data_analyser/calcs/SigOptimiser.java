package com.hartwig.hmftools.data_analyser.calcs;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.calcLinearLeastSquares;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.capValue;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.convertToPercentages;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.copyVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.doubleToStr;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.initVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.lessThan;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sizeToStr;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sumVector;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.data_analyser.types.NmfMatrix;
import com.hartwig.hmftools.data_analyser.types.SampleData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SigOptimiser
{
    // inputs
    private final int mGroupId;
    private List<SampleData> mSamples;
    private List<Integer> mCandiateNewBuckets;
    private List<double[]> mProposedAllocs; // a copy is taken and can be re-adjusted
    private final double[] mStartRatios;

    // initial state
    private int mBucketCount;
    private int mSampleCount;
    private List<Integer> mNonZeroBuckets; // from the non-zero ratios passed in
    private NmfMatrix mSampleCounts; // take from the larger of unallocated + noise and the predefined allocated counts
    private double[] mSampleTotals; // generally each sample's unallocated counts plus noise
    private double[] mProposedTotals; // generally each sample's unallocated counts plus noise
    private double mProposedTotal;
    private NmfMatrix mSampleNoiseRatio; // noise relative to overall sample bucket counts
    private double mCountsTotal; // total of the counts including any applied noise

    // computed state
    private double[] mNoiseRanges;
    private double[] mRatioRanges;
    private double[] mCurrentRatios;
    private double mCurrentAllocPerc;
    private double[] mUnallocatedBucketCounts;
    private double mTotalUnallocCount;
    private List<Integer> mNewBuckets; // added from accepted candidates

    // state
    private boolean mIsValid;
    private boolean mLogVerbose;
    private boolean mHasChanged;

    private static final Logger LOGGER = LogManager.getLogger(SigOptimiser.class);

    public SigOptimiser(int groupId, final List<SampleData> samples, final List<double[]> proposedAllocs,
            final double[] startRatios, List<Integer> candiateNewBuckets)
    {
        mGroupId = groupId;
        mIsValid = true;
        mHasChanged = false;

        if(samples.size() != proposedAllocs.size())
        {
            mIsValid = false;
        }

        mBucketCount = startRatios.length;
        mSampleCount = samples.size();

        mSamples = Lists.newArrayList();
        mSamples.addAll(samples);

        mStartRatios = new double[mBucketCount];
        copyVector(startRatios, mStartRatios);
        mCurrentRatios = new double[mBucketCount];
        copyVector(startRatios, mCurrentRatios);
        mRatioRanges = new double[mBucketCount];
        mNoiseRanges = new double[mBucketCount];

        mNewBuckets = Lists.newArrayList();

        mNonZeroBuckets = Lists.newArrayList();
        for(int b = 0; b < mBucketCount; ++b)
        {
            if(mStartRatios[b] > 0)
                mNonZeroBuckets.add(b);
        }

        mCandiateNewBuckets = Lists.newArrayList();
        mCandiateNewBuckets.addAll(candiateNewBuckets);

        for(Integer bucket : mCandiateNewBuckets)
        {
            if(mNonZeroBuckets.contains(bucket))
            {
                LOGGER.error("candidate bucket({}) part of main set", bucket);
                mIsValid = false;
            }
        }

        mLogVerbose = false;

        // make a matrix from the sample counts
        mSampleCounts = new NmfMatrix(mBucketCount, samples.size());
        mSampleNoiseRatio = new NmfMatrix(mBucketCount, samples.size());

        mProposedAllocs = Lists.newArrayList();
        mProposedAllocs.addAll(proposedAllocs);
        mProposedTotals = new double[mSampleCount];

        mCurrentAllocPerc = 0;
        mUnallocatedBucketCounts = new double[mBucketCount];
        mTotalUnallocCount = 0;

        mCountsTotal = mSampleCounts.sum();
        mSampleTotals = new double[mSampleCount];

        double[][] scData = mSampleCounts.getData();

        for(int s = 0; s < mSampleCount; ++s)
        {
            final SampleData sample = mSamples.get(s);

            final double[] counts = sample.getUnallocBucketCounts();
            final double[] allCounts = sample.getElevatedBucketCounts();
            final double[] noise = sample.getCountRanges();
            final double[] allocNoise = sample.getAllocNoiseCounts();
            final double[] proposedAllocCounts = mProposedAllocs.get(s);

            mProposedTotals[s] = sumVector(proposedAllocCounts);

            for(Integer b : mNonZeroBuckets)
            {
                double maxUnalloc = counts[b] + max(noise[b] - allocNoise[b],0);
                double proposed = proposedAllocCounts[b];

                double sbCount = max(maxUnalloc, proposed);

                if(sbCount == 0)
                {
                    mIsValid = false;
                    LOGGER.error(String.format("grp(%d) sample(%d) invalid bucket(%d) count(elev=%.2f unalloc=%.1f, noise=%.1f/%.1f prop=%.1f)",
                            mGroupId, s, b, allCounts[b], counts[b], noise[b], allocNoise[b], proposedAllocCounts[b]));
                    break;
                }

                scData[b][s] = sbCount;
                mSampleTotals[s] += sbCount;
            }
        }

        mSampleCounts.cacheTranspose();
        mSampleNoiseRatio.cacheTranspose();

        mCountsTotal = mSampleCounts.sum();
        mProposedTotal = sumVector(mProposedTotals);
    }

    public final double[] getFittedRatios() { return mCurrentRatios; }
    public double getAllocPerc() { return mCurrentAllocPerc; }
    public boolean isValid() { return mIsValid; }
    public boolean hasChanged() { return mHasChanged; }
    public List<Integer> getNewBuckets() { return mNewBuckets; }
    public double[] getRevisedSampleAllocCounts(int sampleIndex) { return mProposedAllocs.get(sampleIndex); }
    public double getRevisedSampleAlloc(int sampleIndex) { return mProposedTotals[sampleIndex]; }

    public void setLogVerbose(boolean toggle) { mLogVerbose = toggle; }

    public boolean optimiseBucketRatios()
    {
        if(!mIsValid)
            return false;

        logStats();

        LOGGER.debug(String.format("grp(%d) optimising for samples(%d) buckets(%d) candidates(%d) count(%s unalloc=%s)",
                mGroupId, mSampleCount, mNonZeroBuckets.size(), mCandiateNewBuckets.size(),
                sizeToStr(mCountsTotal), sizeToStr(mTotalUnallocCount)));

        logRatios();

        testCandidateExtraBuckets();

        calcBucketRatioRanges();

        calcOptimalRatios();

        logRatios();

        return true;
    }

    private void calcBucketRatioRanges()
    {
        double[][] nrData = mSampleNoiseRatio.getData();
        double[][] scData = mSampleCounts.getData();

        for(int s = 0; s < mSampleCount; ++s)
        {
            final SampleData sample = mSamples.get(s);

            final double[] allCounts = sample.getElevatedBucketCounts();
            final double[] noise = sample.getCountRanges();

            for(Integer b : mNonZeroBuckets)
            {
                if(allCounts[b] > 0)
                    nrData[b][s] = noise[b]/allCounts[b];
            }
        }

        // calculate a percentage range for noise around each bucket
        final double[] bucketRatios = mCurrentRatios;

        for(int b = 0; b < mBucketCount; ++b)
        {
            if(!mNonZeroBuckets.contains(b))
                continue;

            // for now skip the candidate ones added
            if(mNewBuckets.contains(b))
                continue;

            double percentTotal = 0;
            double noisePercentTotal = 0;
            double countsPercentTotal = 0;

            for(int s = 0; s < mSampleCount; ++s)
            {
                double sampleFactor = mSampleTotals[s] / mCountsTotal; // scale to the significance of this sample to the group
                percentTotal += sampleFactor;

                double noiseRange = nrData[b][s];
                noisePercentTotal += sampleFactor * noiseRange;

                // take any excess sample count as a possible indication of a higher potential ratio
                double countsRatio = scData[b][s] / mProposedTotals[s];

                if(countsRatio > bucketRatios[b])
                {
                    countsPercentTotal += sampleFactor * (countsRatio - bucketRatios[b]);
                }
            }

            double avgNoiseRange = capValue(noisePercentTotal / percentTotal, 0, 0.5);
            double avgCountsRange = capValue(countsPercentTotal / percentTotal, 0, 0.5);

            mRatioRanges[b] = avgCountsRange;
            mNoiseRanges[b] = avgNoiseRange;

            if(mLogVerbose)
            {
                LOGGER.debug(String.format("grp(%d) bucket(%d) ratio(%.4f) range(noise=%.4f counts=%.4f)", mGroupId, b, bucketRatios[b], avgNoiseRange, avgCountsRange));
            }
        }
    }

    private void calcOptimalRatios()
    {
        double[] currentRatios = new double[mBucketCount];
        double[] bestRatios = new double[mBucketCount];
        copyVector(mCurrentRatios, currentRatios);
        copyVector(mCurrentRatios, bestRatios);
        double[] testUnallocCounts = new double[mBucketCount];

        int bucketIncr = 4; // translate this 2*x + 1 test ratios, x negative and x positive

        double bestUnllocTotal = mTotalUnallocCount;

        for(Integer bucket : mNonZeroBuckets)
        {
            double startRatio = bestRatios[bucket]; // start with the existing ratio

            double range = max(mRatioRanges[bucket], mNoiseRanges[bucket]);
            range = min(range, startRatio * 3);

            // double range = mRatioRanges[bucket] / 100; // is a percentage
            double rangeFrag = range / bucketIncr;

            // double bestBucketUnallocTotal = bestUnllocTotal;
            boolean changed = false;
            double bestRatio = 0;

            for(int i = -bucketIncr; i <= bucketIncr; ++i)
            {
                double testRatio = startRatio + i * rangeFrag; // need to use the current normalised ratio as the anchor point

                if(testRatio < 0)
                    continue;

                // always start with the current best set of ratios and add in the trial one
                copyVector(bestRatios, currentRatios);
                currentRatios[bucket] = testRatio;
                convertToPercentages(currentRatios);

                // now test the effect on every sample
                double unallocTotal = calcTotalUnallocatedCounts(currentRatios, testUnallocCounts);

                if (lessThan(unallocTotal, bestUnllocTotal))
                {
                    bestUnllocTotal = unallocTotal;
                    bestRatio = testRatio;
                    changed = true;
                }
            }

            if(changed)
            {
                // keep the best ratios for this bucket and move on
                bestRatios[bucket] = bestRatio;
                convertToPercentages(bestRatios);

                if (mLogVerbose)
                {
                    mHasChanged = true;
                    LOGGER.debug(String.format("grp(%d) bucket(%d) best ratio(%.4f converted=%.4f) vs start(%.4f)",
                            mGroupId, bucket, bestRatio, bestRatios[bucket], mCurrentRatios[bucket]));
                }
            }
        }

        if(bestUnllocTotal < mTotalUnallocCount)
        {
            if(mLogVerbose)
            {
                LOGGER.debug(String.format("grp(%d) lowered unalloc count(%s -> %s)",
                        mGroupId, doubleToStr(mTotalUnallocCount), doubleToStr(bestUnllocTotal)));

                copyVector(bestRatios, mCurrentRatios);
            }

            logStats();
        }
    }

    private void testCandidateExtraBuckets()
    {
        if (mCandiateNewBuckets.isEmpty())
            return;

        /* current logic:
            - some similar groups threw up additional buckets to consider including
            - walk through all samples and extra an implied ratio for each of these buckets
         */
        double[][] scData = mSampleCounts.getData();

        boolean recalcRequired = false;

        for (Integer bucket : mCandiateNewBuckets)
        {
            double[] sampleCounts = new double[mSampleCount];
            double[] sampleRatios = new double[mSampleCount];

            double minSampleRatio = 0;
            double maxSampleRatio = 0;

            for (int s = 0; s < mSampleCount; ++s)
            {
                final SampleData sample = mSamples.get(s);

                double noise = max(sample.getCountRanges()[bucket] - sample.getAllocNoiseCounts()[bucket], 0);
                double sbCount = sample.getUnallocBucketCounts()[bucket] + noise;

                if (sbCount == 0)
                {
                    //if(mLogVerbose)
                    //    LOGGER.debug("grp({}) sample({}) restricts candidate bucket({})", mGroupId, s, bucket);

                    continue;
                }

                // calculate the implied ratio from this count
                double impliedRatio = sbCount / (mProposedTotals[s] + sbCount);

                sampleCounts[s] = sbCount;
                sampleRatios[s] = impliedRatio;

                if (minSampleRatio == 0 || impliedRatio < minSampleRatio)
                {
                    minSampleRatio = impliedRatio;
                }
                else if(impliedRatio > maxSampleRatio)
                {
                    maxSampleRatio = impliedRatio;
                }
            }

            List<Integer> sortedRatioIndices = getSortedVectorIndices(sampleRatios, true);

            // before selecting the optimal ratios, first check that the ratios are somewhat correlated
            // to support the addition of this new bucket to the group
            double lsqRatio = calcLinearLeastSquares(mProposedTotals, sampleCounts);
            double avgRatio = sumVector(sampleCounts) / sumVector(mProposedTotals);
            double medianRatio = 0;

            if(mSampleCount > 10)
            {
                int medianLow = (mSampleCount / 2) - 2;
                int medianHigh = medianLow + 4;
                for(int i = medianLow; i <= medianHigh; ++i)
                {
                    medianRatio += sampleRatios[sortedRatioIndices.get(i)];
                }

                medianRatio /= 5;
            }

            int limitedByCount = 0;
            int limitedByAlloc = 0;
            int limitedByZero = 0;

            // test the median ratio
            double trialRatio = medianRatio;

            double[] adjSampleAllocs = new double[mSampleCount];

            for (int s = 0; s < mSampleCount; ++s)
            {
                if (sampleCounts[s] == 0)
                {
                    adjSampleAllocs[s] = 0;
                    ++limitedByZero;
                    continue;
                }

                double optBucketAlloc = (trialRatio * mProposedTotals[s]) / (1 - trialRatio);

                if(optBucketAlloc < sampleCounts[s])
                {
                    // the ratio would imply all the sample's bucket count could be allocated, so need to limit by the other buckets
                    adjSampleAllocs[s] = mProposedTotals[s] + optBucketAlloc;
                    ++limitedByAlloc;
                }
                else
                {
                    // scale back the proportion of the allocation that could be used
                    double limitFactor = sampleCounts[s] / optBucketAlloc;
                    adjSampleAllocs[s] = (mProposedTotals[s] + sampleCounts[s]) * limitFactor;
                    ++limitedByCount;
                }
            }

            double maxAllocTotal = sumVector(adjSampleAllocs);
            double selectedRatio = trialRatio;

            if(mLogVerbose)
            {
                LOGGER.debug(String.format("grp(%d) candidate bucket(%d) ratios(%.4f min=%.4f max=%.4f lsq=%.4f med=%.4f avg=%.4f) limits(count=%d total=%s zero=%d) with allocTotal(%s) vs total(%s)",
                        mGroupId, bucket, selectedRatio, minSampleRatio, maxSampleRatio, lsqRatio, medianRatio, avgRatio,
                        limitedByCount, limitedByAlloc, limitedByZero, doubleToStr(maxAllocTotal), doubleToStr(mProposedTotal)));
            }

            // test each of the ratios by calculating the allocation gain across each sample to use that ratio
            // any sample with zero counts in this bucket would mean a total loss of its allocation
            // but could still mean an overall gain via the other samples' increases
            int maxLimByCount = 0;
            int maxLimByAlloc = 0;
            int maxLimByNoCount = 0;

            for (Integer ratioIndex : sortedRatioIndices)
            {
                double sampleRatio = sampleRatios[ratioIndex];

                if (sampleRatio == 0)
                    break;

                limitedByCount = 0;
                limitedByAlloc = 0;
                limitedByZero = 0;

                for (int s = 0; s < mSampleCount; ++s)
                {
                    if (sampleCounts[s] == 0)
                    {
                        adjSampleAllocs[s] = 0;
                        ++limitedByZero;
                        continue;
                    }

                    // from solving OptimalBucketCount / Ratio  = TotalExistingAlloc + OptimalBucketCount
                    double optBucketAlloc = (sampleRatio * mProposedTotals[s]) / (1 - sampleRatio);

                    if(optBucketAlloc < sampleCounts[s])
                    {
                        // the ratio would imply all the sample's bucket count could be allocated, so need to limit by the other buckets
                        adjSampleAllocs[s] = mProposedTotals[s] + optBucketAlloc;
                        ++limitedByAlloc;
                    }
                    else
                    {
                        // scale back the proportion of the allocation that could be used
                        double limitFactor = sampleCounts[s] / optBucketAlloc;
                        adjSampleAllocs[s] = (mProposedTotals[s] + sampleCounts[s]) * limitFactor;
                        ++limitedByCount;
                    }
                }

                double allocTotal = sumVector(adjSampleAllocs);

                if (allocTotal > maxAllocTotal)
                {
                    maxAllocTotal = allocTotal;
                    selectedRatio = sampleRatio;
                    maxLimByAlloc = limitedByAlloc;
                    maxLimByCount = limitedByCount;
                    maxLimByNoCount = limitedByZero;
                }
            }

            if(mLogVerbose)
            {
                LOGGER.debug(String.format("grp(%d) candidate bucket(%d) optimal ratio(%.4f min=%.4f max=%.4f) limits(count=%d total=%s zero=%d) with allocTotal(%s) vs total(%s)",
                        mGroupId, bucket, selectedRatio, minSampleRatio, maxSampleRatio,
                        maxLimByCount, maxLimByAlloc, maxLimByNoCount, doubleToStr(maxAllocTotal), doubleToStr(mProposedTotal)));
            }

            if (maxAllocTotal > mProposedTotal)
            {
                LOGGER.debug(String.format("grp(%d) adding candidate bucket(%d) optimal ratio(%.4f) with allocTotal(%s) vs total(%s)",
                        mGroupId, bucket, selectedRatio, doubleToStr(maxAllocTotal), doubleToStr(mProposedTotal)));

                mHasChanged = true;
                recalcRequired = true;
                mNonZeroBuckets.add(bucket);
                mNewBuckets.add(bucket);

                // now add this bucket's sample counts to the totals
                mCurrentRatios[bucket] = selectedRatio;
                convertToPercentages(mCurrentRatios);

                for (int s = 0; s < mSampleCount; ++s)
                {
                    // fill in all missing details for this new bucket
                    double sampleCount = sampleCounts[s];
                    mSampleTotals[s] += sampleCount;
                    mCountsTotal += sampleCount;

                    scData[bucket][s] = sampleCount;

                    // take the lower of the allocation or the actual
                    double optBucketAlloc = (selectedRatio * mProposedTotals[s]) / (1 - selectedRatio);
                    double bucketAlloc = min(sampleCount, optBucketAlloc);
                    mProposedTotals[s] += bucketAlloc;
                    mProposedTotal += bucketAlloc;

                    double[] sampleBucketAllocs = mProposedAllocs.get(s);
                    sampleBucketAllocs[bucket] = bucketAlloc;
                }

                logStats();
            }
        }

        if(recalcRequired)
        {
            convertToPercentages(mCurrentRatios);
        }
    }
    private double calcTotalUnallocatedCounts(final double[] ratios, double[] unallocatedBucketCounts)
    {
        final double[][] scData = mSampleCounts.getData();
        initVector(unallocatedBucketCounts, 0);

        // find the limiting bucket and alloc
        double totalUnallocated = 0;

        for(int s = 0; s < mSampleCount; ++s)
        {
            double minAlloc = 0;

            for(Integer b : mNonZeroBuckets)
            {
                double sampleCount = scData[b][s];
                double alloc = sampleCount / ratios[b];

                if(minAlloc == 0 || alloc < minAlloc)
                    minAlloc = alloc;
            }

            // assert minAlloc > 0

            // translate into counts across the board
            for(Integer b : mNonZeroBuckets)
            {
                double sampleCount = scData[b][s];
                double allocCount = minAlloc * ratios[b];

                double residuals = max(sampleCount - allocCount, 0);
                totalUnallocated += residuals;
                unallocatedBucketCounts[b] += residuals;
            }
        }

        return totalUnallocated;
    }

    private void recalcStats()
    {
        mTotalUnallocCount = calcTotalUnallocatedCounts(mCurrentRatios, mUnallocatedBucketCounts);
    }

    private void logStats()
    {
        recalcStats();

        if(!mLogVerbose)
            return;

        LOGGER.debug(String.format("grp(%d) total(%s) alloc(%s proposed=%s) unalloc(%s perc=%.3f) buckets(%d added=%d)",
                mGroupId, doubleToStr(mCountsTotal), doubleToStr(mCountsTotal - mTotalUnallocCount), doubleToStr(mProposedTotal),
                doubleToStr(mTotalUnallocCount), mTotalUnallocCount/mCountsTotal,
                mNonZeroBuckets.size(), mNewBuckets.size()));
    }

    private void logRatios()
    {
        String bucketRatiosStr = "";
        List<Integer> descBucketRatioIndices = getSortedVectorIndices(mCurrentRatios, false);
        int added = 0;
        for(Integer bucket : descBucketRatioIndices)
        {
            if(mCurrentRatios[bucket] == 0)
                break;

            if (!bucketRatiosStr.isEmpty())
                bucketRatiosStr += ", ";

            bucketRatiosStr += String.format("%d=%.4f", bucket, mCurrentRatios[bucket]);

            ++added;
            if(added >= 20)
                break;
        }

        LOGGER.debug("grp({}) buckets: {}", mGroupId, bucketRatiosStr);

    }

    /*
    private void refineBucketGroupBuckets(BucketGroup bucketGroup)
    {
        if(bucketGroup.getBucketIds().size() <= 5)
            return;

        double minPercentOfMax = 0.001;
        double minRatioTotal = 0.95;
        double bucketRatioTotal = 0;
        double maxBucketRatio = 0;

        double[] bucketRatios = bucketGroup.getBucketRatios();

        List<Integer> sortedRatioIndices = getSortedVectorIndices(bucketRatios, false);
        List<Integer> newBucketsList = Lists.newArrayList();

        for(Integer bucketId : sortedRatioIndices)
        {
            double bucketRatio = bucketRatios[bucketId];

            if(maxBucketRatio == 0)
                maxBucketRatio = bucketRatio;

            if(bucketRatio == 0)
                break;

            if(bucketRatioTotal >= minRatioTotal && bucketRatio < minPercentOfMax * maxBucketRatio)
                break;

            bucketRatioTotal += bucketRatio;
            newBucketsList.add(bucketId);
        }

        if(newBucketsList.size() == bucketGroup.getBucketIds().size())
            return;

        LOGGER.debug(String.format("bg(%d) refined buckets(%d -> %d) ratios(total=%.3f max=%.3f avg=%.3f)",
                bucketGroup.getId(), bucketGroup.getBucketIds().size(), newBucketsList.size(),
                bucketRatioTotal, maxBucketRatio, bucketRatioTotal/newBucketsList.size()));

        bucketGroup.reduceToBucketSet(newBucketsList);
    }

    */
}
