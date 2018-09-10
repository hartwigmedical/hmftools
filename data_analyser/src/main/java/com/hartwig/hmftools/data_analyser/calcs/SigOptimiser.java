package com.hartwig.hmftools.data_analyser.calcs;

import static java.lang.Math.abs;
import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.calcLinearLeastSquares;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.capValue;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.convertToPercentages;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.copyVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.doubleToStr;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.doublesEqual;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.greaterThan;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.initVector;
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
    private List<Integer> mCandidateNewBuckets;
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
    private NmfMatrix mBucketRatioWeights; // weighted frequencies of bucket ratios
    private NmfMatrix mBucketRatioFrequencies; // sample frequencies of bucket ratios
    private double mRefMaxRatio;

    private double mCountsTotal; // total of the counts including any applied noise

    // computed state
    private double[] mNoiseRanges;
    private double[] mRatioRanges;
    private double[] mCurrentRatios;
    private double[] mBucketTotals;
    private double[] mRangesLow;
    private double[] mRangesHigh;
    private double[] mCalcRanges;

    private double mCurrentAllocPerc;
    private double mUnallocTotal;
    private double mAllocTotal;
    private List<Integer> mNewBuckets; // added from accepted candidates
    private List<Integer> mRemovedBuckets; // invalidated from initial buckets
    private List<Integer> mZeroedSamples;

    // state
    private boolean mIsValid;
    private boolean mLogVerbose;
    private boolean mHasChanged;
    private boolean mUseRatioMethod;

    private static double MIN_IMPROVE_PERCENT = 0.005;
    private static int BUCKET_RATIO_SEGMENT_COUNT = 50; // for dividing ratio percents into small segments
    private static double SMALL_RATIO_PERC_CUTOFF = 0.01; // so for a max bucket ratio of 0.36, a minor bucket ratio of 0.0036 would be ignored

    private static final Logger LOGGER = LogManager.getLogger(SigOptimiser.class);

    public SigOptimiser(int groupId, final List<SampleData> samples, final List<double[]> proposedAllocs,
            final double[] startRatios, List<Integer> candiateNewBuckets)
    {
        mGroupId = groupId;
        mIsValid = true;
        mHasChanged = false;

        if (proposedAllocs != null && samples.size() != proposedAllocs.size())
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

        mBucketTotals = new double[mBucketCount];

        mNewBuckets = Lists.newArrayList();
        mRemovedBuckets = Lists.newArrayList();

        mNonZeroBuckets = Lists.newArrayList();
        for (int b = 0; b < mBucketCount; ++b)
        {
            if (mStartRatios[b] > 0)
            {
                mNonZeroBuckets.add(b);
                mRefMaxRatio = max(mRefMaxRatio, mStartRatios[b]);
            }
        }

        mCandidateNewBuckets = Lists.newArrayList();
        mCandidateNewBuckets.addAll(candiateNewBuckets);

        for (Integer bucket : mCandidateNewBuckets)
        {
            if (mNonZeroBuckets.contains(bucket))
            {
                LOGGER.error("candidate bucket({}) part of main set", bucket);
                mIsValid = false;
            }
        }

        mUseRatioMethod = false;
        mLogVerbose = false;

        // make a matrix from the sample counts
        mSampleCounts = new NmfMatrix(mBucketCount, samples.size());
        mSampleNoiseRatio = new NmfMatrix(mBucketCount, samples.size());

        mProposedAllocs = Lists.newArrayList();

        if (proposedAllocs != null)
            mProposedAllocs.addAll(proposedAllocs);

        mProposedTotals = new double[mSampleCount];

        mZeroedSamples = Lists.newArrayList();

        mCurrentAllocPerc = 0;
        mUnallocTotal = 0;
        mAllocTotal = 0;

        // capture the sample counts for all the current (non-zero ratio) buckets
        // the proposed totals are either zero if not yet determined, or whatever the bucket group has currently allocated per sample
        mSampleTotals = new double[mSampleCount];

        double[][] scData = mSampleCounts.getData();

        for (int s = 0; s < mSampleCount; ++s)
        {
            final SampleData sample = mSamples.get(s);

            final double[] counts = sample.getUnallocBucketCounts();
            final double[] allCounts = sample.getElevatedBucketCounts();
            final double[] noise = sample.getCountRanges();
            final double[] allocNoise = sample.getAllocNoiseCounts();

            double[] proposedAllocCounts = null;

            if (proposedAllocs != null)
            {
                // take pre-computed allocations
                proposedAllocCounts = proposedAllocs.get(s);
                mProposedTotals[s] = sumVector(proposedAllocCounts);
            }
            else
            {
                // initialised to zero and will be determined during the optimisation routine
                proposedAllocCounts = new double[mBucketCount];
                mProposedAllocs.add(proposedAllocCounts);
            }

            for (Integer b : mNonZeroBuckets)
            {
                double maxUnalloc = counts[b] + max(noise[b] - allocNoise[b], 0);
                double proposed = proposedAllocCounts[b];

                double sbCount = max(maxUnalloc, proposed);

                if (sbCount == 0)
                    sbCount = max(allCounts[b], noise[b]);

                if (sbCount == 0)
                {

                    mIsValid = false;
                    LOGGER.error(String.format("grp(%d) sample(%d) invalid bucket(%d) count(elev=%.2f unalloc=%.1f, noise=%.1f/%.1f prop=%.1f) propTotal(%.1f)",
                            mGroupId, sample.Id, b, allCounts[b], counts[b], noise[b], allocNoise[b], proposedAllocCounts[b], mProposedTotals[s]));
                    break;
                }

                scData[b][s] = sbCount;
                mSampleTotals[s] += sbCount;
            }
        }

        mSampleCounts.cacheTranspose();
        mSampleNoiseRatio.cacheTranspose();

        mCountsTotal = sumVector(mSampleTotals);
        mProposedTotal = sumVector(mProposedTotals);
    }

    public final double[] getFittedRatios() { return mCurrentRatios; }
    public boolean isValid() { return mIsValid; }
    public boolean hasChanged() { return mHasChanged; }
    public List<Integer> getNewBuckets() { return mNewBuckets; }
    public List<Integer> getNonZeroBuckets() { return mNonZeroBuckets; }
    public double getAllocTotal() { return mAllocTotal; }
    public double[] getRevisedSampleAllocCounts(int sampleIndex) { return mProposedAllocs.get(sampleIndex); }
    public double getRevisedSampleAlloc(int sampleIndex) { return mProposedTotals[sampleIndex]; }
    public double[] getRatioRanges() { return mCalcRanges; }

    public void setLogVerbose(boolean toggle) { mLogVerbose = toggle; }
    public void setUseRatioMethod(boolean toggle) { mUseRatioMethod = toggle; }

    public boolean optimiseBucketRatios()
    {
        if (!mIsValid)
            return false;

        LOGGER.debug(String.format("grp(%d) optimising for samples(%d) buckets(%d) candidates(%d) count(%s alloc=%s unalloc=%s)",
                mGroupId, mSampleCount, mNonZeroBuckets.size(), mCandidateNewBuckets.size(),
                sizeToStr(mCountsTotal), sizeToStr(mAllocTotal), sizeToStr(mUnallocTotal)));

        logStats(false);
        logRatios();

        if(mUseRatioMethod)
        {
            calcBucketDistributions();
        }
        else
        {
            calcBucketRatioRanges();

            testCandidateExtraBuckets();

            if (!mIsValid)
                return false;

            if (!mIsValid)
                return false;

            calcOptimalRatios();
        }

        if (!mIsValid)
            return false;

        if (mHasChanged)
        {
            logStats(false);
            logRatios();
        }

        return true;
    }

    public void calcBucketDistributions()
    {
        mBucketRatioWeights = new NmfMatrix(mBucketCount, BUCKET_RATIO_SEGMENT_COUNT);
        mBucketRatioFrequencies = new NmfMatrix(mBucketCount, BUCKET_RATIO_SEGMENT_COUNT);

        mRangesLow = new double[mBucketCount];
        mRangesHigh = new double[mBucketCount];
        mCalcRanges = new double[mBucketCount];

        List<Integer> bucketIds = Lists.newArrayList();
        bucketIds.addAll(mNonZeroBuckets);
        bucketIds.addAll(mCandidateNewBuckets);

        for(Integer b : bucketIds)
        {
            for(final SampleData sample : mSamples)
                mBucketTotals[b] += sample.getElevatedBucketCounts()[b];
        }

        double[] ratioIncrements = new double[mBucketCount];
        initVector(ratioIncrements, 1 / (double) BUCKET_RATIO_SEGMENT_COUNT);

        assignBucketRatioFrequencies(bucketIds, ratioIncrements, mBucketTotals);

        // calculate a cache a mean and max for each bucket from the weighted sample count ratios
        double[] meanBucketRatios = new double[mBucketCount];
        double[] maxBucketRatios = new double[mBucketCount];

        calcBucketDistributionStats(bucketIds, ratioIncrements, meanBucketRatios, maxBucketRatios);

        // scale the frequencies to a more appropriate set of incremental ratio values
        for(Integer b : bucketIds)
        {
            if (meanBucketRatios[b] == 0)
                continue;

            double refRatio = max(meanBucketRatios[b], maxBucketRatios[b]);

            if(refRatio >= 0.5)
                ratioIncrements[b] = 1 / (double) BUCKET_RATIO_SEGMENT_COUNT;
            else if(refRatio >= 0.25)
                ratioIncrements[b] = 0.5 / (double) BUCKET_RATIO_SEGMENT_COUNT;
            else if(refRatio >= 0.1)
                ratioIncrements[b] = 0.2 / (double) BUCKET_RATIO_SEGMENT_COUNT;
            else if(refRatio >= 0.05)
                ratioIncrements[b] = 0.1 / (double) BUCKET_RATIO_SEGMENT_COUNT;
            else
                ratioIncrements[b] = 0.05 / (double) BUCKET_RATIO_SEGMENT_COUNT;
        }

        assignBucketRatioFrequencies(bucketIds, ratioIncrements, mBucketTotals);

        calcBucketDistributionStats(bucketIds, ratioIncrements, meanBucketRatios, maxBucketRatios);

        // assess the results
        double maxMeanBucketRatio = 0;

        for(Integer b : bucketIds)
        {
            if(meanBucketRatios[b] > maxMeanBucketRatio)
                maxMeanBucketRatio = meanBucketRatios[b];
        }

        for(Integer b : bucketIds)
        {
            examineBucketDistribution(b, meanBucketRatios[b], ratioIncrements[b]);
        }

        // finally update the final ratios
        mCountsTotal = sumVector(mSampleTotals);
        mProposedTotal = sumVector(mProposedTotals);

        convertToPercentages(mCurrentRatios);

        logStats(true);
    }

    private void assignBucketRatioFrequencies(final List<Integer> bucketIds, final double[] ratioIncrements, final double[] bucketTotals)
    {
        double[][] brwData = mBucketRatioWeights.getData();
        double[][] brfData = mBucketRatioFrequencies.getData();

        // re-init data
        for(Integer b : bucketIds)
        {
            for (int i = 0; i < BUCKET_RATIO_SEGMENT_COUNT; ++i)
            {
                brwData[b][i] = 0;
                brfData[b][i] = 0;
            }
        }

        for(final SampleData sample : mSamples)
        {
            double sampleTotal = sample.getTotalCount();

            final double[] counts = sample.getUnallocBucketCounts();
            // final double[] sampleCounts = sample.getElevatedBucketCounts();
            final double[] noise = sample.getCountRanges();
            final double[] allocNoise = sample.getAllocNoiseCounts();

            double sampleBgTotal = 0; // total count across all buckets just in this sig
            for(Integer b : bucketIds)
            {
                double sbCount = counts[b];

                if(sbCount == 0)
                    sbCount = max(noise[b] - allocNoise[b], 0);

                sampleBgTotal += sbCount;
            }

            for(Integer b : bucketIds)
            {
                double sbCount = counts[b];

                if(sbCount == 0)
                    sbCount = max(noise[b] - allocNoise[b], 0);

                double rawBucketRatio = sbCount / sampleBgTotal;
                int ratioIndex = (int)round(rawBucketRatio / ratioIncrements[b]);

                if(ratioIndex < 0 || ratioIndex >= BUCKET_RATIO_SEGMENT_COUNT)
                {
                    LOGGER.debug(String.format("sample(%d) bucket(%d) skipping ratioIndex(%d) outside range from sbCount(%s) sampleTotal(%s) rawRatio(%.4f)",
                            sample.Id, b, ratioIndex, sizeToStr(sbCount), sizeToStr(sampleTotal), rawBucketRatio));
                    continue;
                }

                double sampleWeight = sbCount / bucketTotals[b];
                brwData[b][ratioIndex] += sampleWeight;
                brfData[b][ratioIndex] += 1;
            }
        }
    }

    private void calcBucketDistributionStats(final List<Integer> bucketIds, final double[] ratioIncrements, double[] meanBucketRatios, double[] maxBucketRatios)
    {
        double[][] brwData = mBucketRatioWeights.getData();

        for(Integer b : bucketIds)
        {
            double bucketRatioTotal = 0;
            double weightTotal = 0;
            double maxBucketWeight = 0;

            for(int i = 0; i < BUCKET_RATIO_SEGMENT_COUNT; ++i)
            {
                double sampleWeight = brwData[b][i];
                double bucketRatio = i * ratioIncrements[b];
                bucketRatioTotal += bucketRatio * sampleWeight;
                weightTotal += sampleWeight;

                // record the highest weighted ratio
                if (sampleWeight > maxBucketWeight)
                    maxBucketWeight = sampleWeight;

                // record the highest ratio of any significance
                if(maxBucketWeight > 0 && sampleWeight >= 0.1 * maxBucketWeight)
                    maxBucketRatios[b] = bucketRatio;
            }

            if(weightTotal == 0)
                continue;

            meanBucketRatios[b] = bucketRatioTotal / weightTotal;

            // LOGGER.debug(String.format("grp(%d) bucket(%d) mean(%.3f) max(%.3f) ratioInc(%.3f)",
            //         mGroupId, b, meanBucketRatios[b], maxBucketRatios[b], ratioIncrements[b]));
        }
    }

    private static double BUCKET_RATIO_RANGE_PERCENTILE = 0.8;
    private static double BUCKET_RANGE_REQ_SAMPLE_PERC = 0.75;
    private static double BUCKET_RANGE_SCORE_THRESHOLD = 0.05; // 0.5 allocation of samples and >= 10% of max ratio

    private void examineBucketDistribution(int bucket, double origMeanRatio, double ratioIncrement)
    {
        double lastSampleWeight = 0;

        double[][] brwData = mBucketRatioWeights.getData();
        double[][] brfData = mBucketRatioFrequencies.getData();

        double maxBucketWeight = 0;
        double cumulativeFreq = 0;

        double freqTotal = sumVector(mBucketRatioFrequencies.getRow(bucket)); // should be number of samples

        double startPercentile = (1 - BUCKET_RATIO_RANGE_PERCENTILE) * 0.5;
        double startFreqTotal = startPercentile * freqTotal;
        double endFreqTotal = (1 - startPercentile) * freqTotal;

        // now search for a peak and determine a range around it
        int distStartIndex = -1;
        int distEndIndex = -1;

        // find the percentile range based on sample frequencies
        for(int i = 0; i < BUCKET_RATIO_SEGMENT_COUNT; ++i)
        {
            double sampleWeight = brwData[bucket][i];

            if (sampleWeight > maxBucketWeight)
                maxBucketWeight = sampleWeight;

            cumulativeFreq += brfData[bucket][i];

            if(distStartIndex == -1 && cumulativeFreq >= startFreqTotal)
                distStartIndex = i;
            else if(distEndIndex == -1 && cumulativeFreq >= endFreqTotal)
                distEndIndex = i;
        }

        if(distStartIndex == -1 || distEndIndex == -1 || distStartIndex >= distEndIndex)
        {
            LOGGER.error("grp({}) bucket({}) invalid range indices", mGroupId, bucket, distStartIndex, distEndIndex);
            return;
        }

        int maxIndexRange = 10; // make a proportion of the segement count or relate to % of bucket ratio
        int distRange = distEndIndex - distStartIndex;
        if(distRange > maxIndexRange)
        {
            int excess = (int)floor((distRange - maxIndexRange) * 0.5);
            distStartIndex += excess;
            distEndIndex -= excess;
        }

        double bucketRatioTotal = 0;
        double weightTotal = 0;
        // double invalidChangeTotal = 0;
        double rangeBoundMax = 0;
        double maxWeight = 0;
        // int maxWeightIndex = 0;
        for(int i = distStartIndex; i <= distEndIndex; ++i)
        {
            double sampleWeight = brwData[bucket][i];
            double bucketRatio = i * ratioIncrement;
            bucketRatioTotal += bucketRatio * sampleWeight;
            weightTotal += sampleWeight;

            if(sampleWeight > maxWeight)
            {
                maxWeight = sampleWeight;
                rangeBoundMax = bucketRatio;
            }

            /*
            if(i == 0)
            {
                lastSampleWeight = sampleWeight;
                continue;
            }

            if(lastSampleWeight == 0 && sampleWeight == 0)
                continue;

            double swChange = sampleWeight - lastSampleWeight;
            boolean isSignificantWeight = sampleWeight > 0.25 * maxBucketWeight;

            if(swChange > 0 && i > maxWeightIndex && isSignificantWeight && distEndIndex == -1)
            {
                // heading up after the peak within percentile range, may point to additional peaks
                invalidChangeTotal += swChange;
            }
            else if(swChange < 0 && i < maxWeightIndex && isSignificantWeight && distStartIndex >= 0)
            {
                // heading down before the peak within percentile range, may point to additional peaks
                invalidChangeTotal += abs(swChange);
            }

            lastSampleWeight = sampleWeight;
            */
        }

        double rangeBoundMean = bucketRatioTotal / weightTotal;
        double startRatio = distStartIndex >= 0 ? distStartIndex * ratioIncrement : 0;
        double endRatio = distEndIndex >= 0 ? distEndIndex * ratioIncrement : 0;

        mRangesLow[bucket] = startRatio;
        mRangesHigh[bucket] = endRatio;
        mCurrentRatios[bucket] = rangeBoundMean;

        double minRange = min(rangeBoundMean - startRatio, endRatio - rangeBoundMean);

        double allocScore = (rangeBoundMean / mRefMaxRatio) * weightTotal;
        boolean isCandidate = mCandidateNewBuckets.contains(bucket);

        if(mLogVerbose)
        {
            double countsTotal = weightTotal * mBucketTotals[bucket];

            LOGGER.debug(String.format("grp(%d) bucket(%d %s) score(%.3f) mean(%.3f raw=%.3f) range(%.3f -> %.3f min=%.3f indx=%d -> %d) max(%.3f) total(%.3f %s of %s)",
                    mGroupId, bucket, isCandidate ? "cand" : "init", allocScore, rangeBoundMean, origMeanRatio, startRatio, endRatio, minRange,
                    distStartIndex, distEndIndex, rangeBoundMax, weightTotal, sizeToStr(countsTotal), sizeToStr(mBucketTotals[bucket])));
        }

        if(!isCandidate)
        {
            // always take the range-bound calculate mean
            mCurrentRatios[bucket] = rangeBoundMean;
            mCalcRanges[bucket] = minRange;

            mHasChanged = true;
            return;
        }

        // include buckets if a sufficient proportion of sample counts are included within the ranges
        // and which aren't too small relative to the main bucket(s) to be examined with any precision
        boolean bucketValid = distStartIndex > 0 && allocScore >= BUCKET_RANGE_SCORE_THRESHOLD;

        if(!bucketValid)
            return;

        mNewBuckets.add(bucket);
        mNonZeroBuckets.add(bucket);
        mHasChanged = true;

        double[][] scData = mSampleCounts.getData();

        // now add this bucket's sample counts to the totals
        mCurrentRatios[bucket] = rangeBoundMean;
        mCalcRanges[bucket] = minRange;

        for (int s = 0; s < mSampleCount; ++s)
        {
            // fill in all missing details for this new bucket
            final SampleData sample = mSamples.get(s);
            double sampleCount = sample.getElevatedBucketCounts()[bucket];
            mSampleTotals[s] += sampleCount;
            scData[bucket][s] = sampleCount;

            // take the lower of the allocation or the actual
            double optBucketAlloc = (rangeBoundMean * mProposedTotals[s]) / (1 - rangeBoundMean);
            double bucketAlloc = min(sampleCount, optBucketAlloc);
            mProposedTotals[s] += bucketAlloc;

            double[] sampleBucketAllocs = mProposedAllocs.get(s);
            sampleBucketAllocs[bucket] = bucketAlloc;
        }
    }

    private void calcBucketRatioRanges()
    {
        double[][] nrData = mSampleNoiseRatio.getData();
        double[][] scData = mSampleCounts.getData();

        for (int s = 0; s < mSampleCount; ++s)
        {
            final SampleData sample = mSamples.get(s);

            final double[] allCounts = sample.getElevatedBucketCounts();
            final double[] noise = sample.getCountRanges();

            for (Integer b : mNonZeroBuckets)
            {
                if (allCounts[b] > 0)
                    nrData[b][s] = noise[b] / allCounts[b];
            }
        }

        // calculate a percentage range for noise around each bucket
        final double[] bucketRatios = mCurrentRatios;

        for (int b = 0; b < mBucketCount; ++b)
        {
            if (!mNonZeroBuckets.contains(b))
                continue;

            // for now skip the candidate ones added
            if (mNewBuckets.contains(b))
                continue;

            double percentTotal = 0;
            double noisePercentTotal = 0;
            double countsPercentTotal = 0;

            for (int s = 0; s < mSampleCount; ++s)
            {
                double sampleFactor = mSampleTotals[s] / mCountsTotal; // scale to the significance of this sample to the group
                percentTotal += sampleFactor;

                double noiseRange = nrData[b][s];
                noisePercentTotal += sampleFactor * noiseRange;

                // take any excess sample count as a possible indication of a higher potential ratio
                double countsRatio = scData[b][s] / mProposedTotals[s];

                if (countsRatio > bucketRatios[b])
                {
                    countsPercentTotal += sampleFactor * (countsRatio - bucketRatios[b]);
                }
            }

            double avgNoiseRange = capValue(noisePercentTotal / percentTotal, 0, 0.25);
            double avgCountsRange = capValue(countsPercentTotal / percentTotal, 0, 0.5);

            mRatioRanges[b] = avgCountsRange;
            mNoiseRanges[b] = avgNoiseRange;

            if (mLogVerbose)
            {
                LOGGER.debug(String.format("grp(%d) bucket(%d) ratio(%.4f) range(noise=%.4f counts=%.4f)",
                        mGroupId, b, bucketRatios[b], avgNoiseRange, avgCountsRange));
            }
        }
    }

    private void calcOptimalRatios()
    {
        double[] currentRatios = new double[mBucketCount];
        double[] bestRatios = new double[mBucketCount];
        copyVector(mCurrentRatios, currentRatios);
        copyVector(mCurrentRatios, bestRatios);
        double[] testSampleAllocs = new double[mSampleCount];

        int bucketIncr = 5;

        double bestAllocTotal = 0;
        List<Integer> bucketsToRemove = Lists.newArrayList();
        List<Double> testRatios = Lists.newArrayList();

        for (Integer bucket : mNonZeroBuckets)
        {
            double startRatio = bestRatios[bucket]; // start with the existing ratio

            double range = max(mRatioRanges[bucket], mNoiseRanges[bucket]);

            boolean changed = false;
            double bestRatio = 0;
            double testMin = 1;
            double testMax = 0;
            double bucketBestAllocTotal = 0;
            double bucketBestRatio = 0;

            // test a range of values around the starting ratio
            // and repeat the process with a tighter range centered around the new best ration the second time through
            for (int j = 0; j < 2; ++j)
            {
                testRatios.clear();

                if (j == 1)
                {
                    if (bucketBestRatio == startRatio)
                    {
                        // try again but with a much tighter increment
                        range *= 0.2;
                    }
                    else
                    {
                        // refine the search
                        range = abs(startRatio - bucketBestRatio);
                        startRatio = bucketBestRatio;
                    }
                }

                if (startRatio - range <= 0)
                    testRatios.add(0.0);

                double rangeFrag;

                if (range == 0)
                    rangeFrag = startRatio * (0.5 / (double) bucketIncr); // eg 10% for 5 increments
                else
                    rangeFrag = min(range, startRatio) / bucketIncr;

                for (int i = -bucketIncr; i <= bucketIncr + 1; ++i)
                {
                    double testRatio = startRatio + i * rangeFrag;

                    if (testRatio > 0)
                        testRatios.add(testRatio);
                }

                for (Double testRatio : testRatios)
                {
                    testMax = max(testMax, testRatio);
                    testMin = min(testMin, testRatio);

                    // always start with the current best set of ratios and add in the trial one
                    copyVector(bestRatios, currentRatios);
                    currentRatios[bucket] = testRatio;
                    convertToPercentages(currentRatios);

                    // now test the effect on every sample
                    // double unallocTotal = calcTotalUnallocatedCounts(currentRatios, testUnallocCounts, null);
                    double allocTotal = calcTotalAllocatedCounts(currentRatios, testSampleAllocs, null);

                    if (greaterThan(allocTotal, bucketBestAllocTotal))
                    {
                        bucketBestAllocTotal = allocTotal;
                        bucketBestRatio = testRatio;
                    }

                    if (greaterThan(allocTotal, bestAllocTotal))
                    {
                        bestAllocTotal = allocTotal;
                        bestRatio = testRatio;
                        changed = true;
                    }
                }
            }

            if (changed)
            {
                // keep the best ratios for this bucket and move on
                bestRatios[bucket] = bestRatio;
                convertToPercentages(bestRatios);

                mHasChanged = true;
                LOGGER.debug(String.format("grp(%d) bucket(%d) best ratio(%.4f converted=%.4f) vs start(%.4f) testRange(%.4f -> %.4f)",
                        mGroupId, bucket, bestRatio, bestRatios[bucket], mCurrentRatios[bucket], testMin, testMax));

                if (bestRatio == 0)
                    bucketsToRemove.add(bucket);
            }
        }

        for (Integer bucket : bucketsToRemove)
        {
            mNonZeroBuckets.remove(bucket);
            mRemovedBuckets.add(bucket);

            for (int s = 0; s < mSampleCount; ++s)
            {
                double[] sampleBucketAllocs = mProposedAllocs.get(s);
                sampleBucketAllocs[bucket] = 0;
            }
        }

        if (bestAllocTotal > mAllocTotal)
        {
            LOGGER.debug(String.format("grp(%d) ratio optimisation increased alloc count(%s -> %s)",
                    mGroupId, doubleToStr(mAllocTotal), doubleToStr(bestAllocTotal)));

            copyVector(bestRatios, mCurrentRatios);

            logStats(true);
        }
    }

    private void testCandidateExtraBuckets()
    {
        if (mCandidateNewBuckets.isEmpty())
            return;

        /* current logic:
            - some similar groups threw up additional buckets to consider including
            - walk through all samples and extract an implied ratio for each of these buckets
         */
        double[][] scData = mSampleCounts.getData();

        boolean recalcRequired = false;
        double minRatioRangePerc = 0.2;

        for (Integer bucket : mCandidateNewBuckets)
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
                else if (impliedRatio > maxSampleRatio)
                {
                    maxSampleRatio = impliedRatio;
                }
            }

            // before selecting the optimal ratios, first check that the ratios are somewhat correlated
            // to support the addition of this new bucket to the group
            double lsqRatio = calcLinearLeastSquares(mProposedTotals, sampleCounts);
            double avgRatio = sumVector(sampleCounts) / sumVector(mProposedTotals);

            if (mSampleCount < 10)
            {
                // take the bucket only if the correlation is very high
                // in which case the bucket should have already been included
                continue;
            }

            // new test: 80% of samples are within X% of each other
            List<Integer> sortedRatioIndices = getSortedVectorIndices(sampleRatios, true);

            double percentRange = 0.8;
            double maxRatioDiff = 0.2; // eg 0.2 to 0.16 or 0.24
            double percLow = ((1 - percentRange) * 0.5);
            double percHigh = 1 - percLow;
            int indexLowPerc = (int) round(mSampleCount * percLow);
            int indexHighPerc = (int) round(mSampleCount * percHigh);
            double ratioLowPerc = sampleRatios[sortedRatioIndices.get(indexLowPerc)];
            double ratioHighPerc = sampleRatios[sortedRatioIndices.get(indexHighPerc)];

            int medianIndex = (int) round(mSampleCount * 0.5);
            double medianRatio = sampleRatios[sortedRatioIndices.get(medianIndex)];

            if (medianRatio == 0)
                continue;

            if ((medianRatio - ratioLowPerc) / medianRatio > maxRatioDiff)
                continue;

            if ((ratioHighPerc - medianRatio) / medianRatio > maxRatioDiff)
                continue;

            double trialRatio = medianRatio;

            int limitedByCount = 0;
            int limitedByAlloc = 0;
            int limitedByZero = 0;

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

                if (optBucketAlloc < sampleCounts[s])
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

            if (mLogVerbose)
            {
                LOGGER.debug(String.format("grp(%d) candidate bucket(%d) ratios(%.4f min=%.4f max=%.4f lsq=%.4f med=%.4f avg=%.4f lowp=%.4f highp=%.4f) limits(count=%d total=%s zero=%d) with allocTotal(%s) vs total(%s)",
                        mGroupId, bucket, selectedRatio, minSampleRatio, maxSampleRatio, lsqRatio, medianRatio, avgRatio, ratioLowPerc, ratioHighPerc,
                        limitedByCount, limitedByAlloc, limitedByZero, doubleToStr(maxAllocTotal), doubleToStr(mProposedTotal)));
            }

            /*
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
            */

            double percIncrease = (maxAllocTotal - mAllocTotal) / mAllocTotal;

            if (percIncrease >= MIN_IMPROVE_PERCENT)
            {
                LOGGER.debug(String.format("grp(%d) adding candidate bucket(%d) optimal ratio(%.4f) with alloc improve(%.3f %s prev=%s) ratios(min=%.4f max=%.4f lsq=%.4f med=%.4f avg=%.4f)",
                        mGroupId, bucket, selectedRatio, percIncrease, doubleToStr(maxAllocTotal), doubleToStr(mAllocTotal),
                        minSampleRatio, maxSampleRatio, lsqRatio, trialRatio, avgRatio));

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
                    scData[bucket][s] = sampleCount;

                    // take the lower of the allocation or the actual
                    double optBucketAlloc = (selectedRatio * mProposedTotals[s]) / (1 - selectedRatio);
                    double bucketAlloc = min(sampleCount, optBucketAlloc);
                    mProposedTotals[s] += bucketAlloc;

                    double[] sampleBucketAllocs = mProposedAllocs.get(s);
                    sampleBucketAllocs[bucket] = bucketAlloc;
                }

                mCountsTotal = sumVector(mSampleTotals);
                mProposedTotal = sumVector(mProposedTotals);

                logStats(true);
            }
        }

        if (recalcRequired)
        {
            convertToPercentages(mCurrentRatios);
        }
    }

    private double calcTotalAllocatedCounts(final double[] ratios, double[] allocatedSampleCounts, List<Integer> zeroedSamples)
    {
        // goes through each sample and applies the input ratios to calculate what can be allocated
        // return data:
        // - allocated bucket counts - per-bucket counts with the input ratios
        // - zeroed samples - any sample which cannot allocate anthing with the input ratios
        // - total allocated count across all samples

        final double[][] scData = mSampleCounts.getData();
        initVector(allocatedSampleCounts, 0);

        // find the limiting bucket and alloc
        double totalAllocated = 0;

        for (int s = 0; s < mSampleCount; ++s)
        {
            double minAlloc = 0;

            for (Integer b : mNonZeroBuckets)
            {
                double sampleCount = scData[b][s];
                double alloc = sampleCount / ratios[b];

                if (minAlloc == 0 || alloc < minAlloc)
                    minAlloc = alloc;
            }

            if (minAlloc == 0 && zeroedSamples != null)
                zeroedSamples.add(s);

            allocatedSampleCounts[s] = minAlloc;

            // translate into counts across the board
            for (Integer b : mNonZeroBuckets)
            {
                double allocCount = minAlloc * ratios[b];
                totalAllocated += allocCount;
            }
        }

        return totalAllocated;
    }

    private void recalcStats()
    {
        mZeroedSamples.clear();

        double ratioSum = sumVector(mCurrentRatios);

        if (!doublesEqual(ratioSum, 1))
        {
            LOGGER.debug(String.format("grp(%d) invalid ratioSum(%.6f", mGroupId, ratioSum));
            logRatios();
            mIsValid = false;
            return;
        }

        mAllocTotal = calcTotalAllocatedCounts(mCurrentRatios, mProposedTotals, mZeroedSamples);

        if (greaterThan(mAllocTotal, mCountsTotal))
        {
            LOGGER.debug(String.format("grp(%d) allocTotal(%.3f) exceeds sampleCountsTotal(%.3f)", mGroupId, mAllocTotal, mCountsTotal));
            mIsValid = false;
            return;
        }

        mUnallocTotal = mCountsTotal - mAllocTotal;
    }

    private void logStats(boolean checkVerbose)
    {
        recalcStats();

        if (!mLogVerbose && !checkVerbose)
            return;

        LOGGER.debug(String.format("grp(%d) total(%s) alloc(%s perc=%.3f proposed=%s) unalloc(%s perc=%.3f) buckets(%d add=%d remove=%d) samples(%d zeroed=%d)",
                mGroupId, doubleToStr(mCountsTotal), doubleToStr(mAllocTotal), mAllocTotal / mCountsTotal, doubleToStr(mProposedTotal),
                doubleToStr(mUnallocTotal), mUnallocTotal / mCountsTotal,
                mNonZeroBuckets.size(), mNewBuckets.size(), mRemovedBuckets.size(), mSampleCount, mZeroedSamples.size()));
    }

    private void logRatios()
    {
        String bucketRatiosStr = "";
        List<Integer> descBucketRatioIndices = getSortedVectorIndices(mCurrentRatios, false);
        int added = 0;
        for (Integer bucket : descBucketRatioIndices)
        {
            if (mCurrentRatios[bucket] == 0)
                break;

            if (!bucketRatiosStr.isEmpty())
                bucketRatiosStr += ", ";

            bucketRatiosStr += String.format("%d=%.4f", bucket, mCurrentRatios[bucket]);

            ++added;
            if (added >= 20)
                break;
        }

        LOGGER.debug("grp({}) buckets: {}", mGroupId, bucketRatiosStr);
    }


}
