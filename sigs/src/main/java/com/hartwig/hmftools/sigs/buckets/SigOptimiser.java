package com.hartwig.hmftools.sigs.buckets;

import static java.lang.Math.abs;
import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.sigs.DataUtils.doubleToStr;
import static com.hartwig.hmftools.common.sigs.DataUtils.doublesEqual;
import static com.hartwig.hmftools.common.sigs.DataUtils.greaterThan;
import static com.hartwig.hmftools.common.sigs.DataUtils.sizeToStr;
import static com.hartwig.hmftools.common.sigs.SigUtils.convertToPercentages;
import static com.hartwig.hmftools.common.utils.VectorUtils.copyVector;
import static com.hartwig.hmftools.common.utils.VectorUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.common.utils.VectorUtils.initVector;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;

import java.util.HashMap;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Matrix;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SigOptimiser
{
    // inputs
    private final int mGroupId;
    private List<SampleData> mAllSamples;
    private List<Integer> mCandidateNewBuckets;
    private List<double[]> mProposedAllocs; // a copy is taken and can be re-adjusted
    private final double[] mStartRatios;
    private boolean mAssessmentOnly;

    private double mRatioRangePercentile; // percentage of middle sample's ratio data to use

    // initial state
    private int mBucketCount;
    private int mAllSampleCount;
    private List<Integer> mInitialBuckets; // from the non-zero ratios passed in
    private Matrix mSampleCounts; // take from the larger of unallocated + noise and the predefined allocated counts
    private double[] mSampleTotals; // generally each sample's unallocated counts plus noise
    private double[] mProposedTotals; // generally each sample's unallocated counts plus noise
    private double mProposedTotal;
    private double mRefMaxRatio;

    private double mCountsTotal; // total of the counts including any applied noise

    // computed values
    private List<SampleData> mSamples; // only those which are primarily covered by this sig and so considered pure for optimisation
    private int mSampleCount;
    private double[] mCurrentRatios;
    private double[] mRangeMeanRatios; // weighted mean within the distribution range
    private double[] mRangesLow;
    private double[] mRangesHigh;
    private double[] mCalcRanges;
    private double[] mUnallocBucketTotals; // from pure samples' unallocated counts
    private Matrix mSampleNoiseRatio; // noise relative to overall sample bucket counts
    private Matrix mBucketRatioWeights; // weighted frequencies of bucket ratios
    private Matrix mBucketRatioFrequencies; // sample frequencies of bucket ratios

    private Matrix mSampleSigContrib; // sample sig contributions
    private Matrix mSampleRatios; // sample contributions turned into ratios

    private double mUnallocTotal;
    private double mAllocTotal;
    private List<Integer> mActiveBuckets; // buckets which will form part of the signature at the end
    private List<Integer> mNewBuckets; // added from accepted candidates
    private List<Integer> mRemovedBuckets; // invalidated from initial buckets
    private List<Integer> mZeroedSamples;

    private boolean mCacheBucketInfo;
    private HashMap<Integer,String> mBucketInfo; // for external logging

    // state
    private boolean mIsValid;
    private boolean mLogVerbose;
    private boolean mHasChanged;

    private static final int MIN_PURE_SAMPLE_COUNT = 5;
    public static final double SMALL_RATIO_PERC_CUTOFF = 0.01; // so for a max bucket ratio of 0.36, a minor bucket ratio of 0.0036 would be ignored
    private static final double DEFAULT_RATIO_RANGE_PERCENTILE = 0.6;
    public static final double BUCKET_RANGE_MAX_PERCENT = 0.20; // max that a range can be as a percentage of its ratio, works up and down
    public static final double SAMPLE_PURE_SIG_PERCENT = 0.8; // used to decide whether a sample can aid in optimising sig ratios and ranges

    private static final Logger LOGGER = LogManager.getLogger(SigOptimiser.class);

    public SigOptimiser(int groupId, final List<SampleData> samples, final List<double[]> proposedAllocs,
            final double[] startRatios, List<Integer> candidateNewBuckets)
    {
        mGroupId = groupId;
        mIsValid = true;
        mHasChanged = false;
        mAssessmentOnly = false;

        mRatioRangePercentile = DEFAULT_RATIO_RANGE_PERCENTILE;

        if (proposedAllocs != null && samples.size() != proposedAllocs.size())
        {
            mIsValid = false;
        }

        mBucketCount = startRatios.length;
        mAllSampleCount = samples.size();

        mAllSamples = Lists.newArrayList();
        mAllSamples.addAll(samples);

        mBucketInfo = Maps.newHashMap();

        mStartRatios = new double[mBucketCount];
        copyVector(startRatios, mStartRatios);
        mCurrentRatios = new double[mBucketCount];
        copyVector(startRatios, mCurrentRatios);
        mRangeMeanRatios = new double[mBucketCount];

        mNewBuckets = Lists.newArrayList();
        mRemovedBuckets = Lists.newArrayList();
        mInitialBuckets = Lists.newArrayList();
        mActiveBuckets = Lists.newArrayList();

        for (int b = 0; b < mBucketCount; ++b)
        {
            if (mStartRatios[b] > 0)
            {
                mInitialBuckets.add(b);
                mRefMaxRatio = max(mRefMaxRatio, mStartRatios[b]);
            }
        }

        mActiveBuckets.addAll(mInitialBuckets);

        mCandidateNewBuckets = Lists.newArrayList();
        mCandidateNewBuckets.addAll(candidateNewBuckets);

        for (Integer bucket : mCandidateNewBuckets)
        {
            if (mActiveBuckets.contains(bucket))
            {
                LOGGER.error("candidate bucket({}) part of main set", bucket);
                mIsValid = false;
            }
        }

        mSamples = Lists.newArrayList();

        findPureSamples();

        mLogVerbose = false;
        mCacheBucketInfo = false;

        // make a matrix from the sample counts
        mSampleCounts = new Matrix(mBucketCount, mAllSampleCount);
        mSampleNoiseRatio = new Matrix(mBucketCount, mAllSampleCount);

        mProposedAllocs = Lists.newArrayList();

        if (proposedAllocs != null)
            mProposedAllocs.addAll(proposedAllocs);

        mProposedTotals = new double[mAllSampleCount];

        mZeroedSamples = Lists.newArrayList();

        mUnallocTotal = 0;
        mAllocTotal = 0;

        // capture the sample counts for all the current (non-zero ratio) buckets
        // the proposed totals are either zero if not yet determined, or whatever the bucket group has currently allocated per sample
        mSampleTotals = new double[mAllSampleCount];

        double[][] scData = mSampleCounts.getData();

        for (int s = 0; s < mAllSampleCount; ++s)
        {
            final SampleData sample = mAllSamples.get(s);

            final double[] counts = sample.getUnallocBucketCounts();
            final double[] allCounts = sample.getElevatedBucketCounts();
            final double[] noise = sample.getNoiseCounts();
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

            for (Integer b : mActiveBuckets)
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

        mUnallocBucketTotals = new double[mBucketCount];

        List<Integer> allBuckets = Lists.newArrayList();
        allBuckets.addAll(mActiveBuckets);
        allBuckets.addAll(mCandidateNewBuckets);

        for(final SampleData sample : mSamples)
        {
            final double[] counts = sample.getUnallocBucketCounts();

            for (Integer b : allBuckets)
            {
                mUnallocBucketTotals[b] += counts[b];
            }
        }

        mSampleCounts.cacheTranspose();
        mSampleNoiseRatio.cacheTranspose();

        mCountsTotal = sumVector(mSampleTotals);
        mProposedTotal = sumVector(mProposedTotals);

        mRangesLow = new double[mBucketCount];
        mRangesHigh = new double[mBucketCount];
        mCalcRanges = new double[mBucketCount];
    }

    private void findPureSamples()
    {
        for(final SampleData sample : mAllSamples)
        {
            final double[] allCounts = sample.getElevatedBucketCounts();
            double sampleSigTotal = 0;

            for (Integer b : mActiveBuckets)
            {
                sampleSigTotal += allCounts[b];
            }

            if(sampleSigTotal / sample.getElevatedCount() > SAMPLE_PURE_SIG_PERCENT)
            {
                mSamples.add(sample);
            }
        }

        mSampleCount = mSamples.size();

        if(mSampleCount == 0)
        {
            LOGGER.debug("grp({}) has no pure samples", mGroupId);
            mIsValid = false;
        }
        else
        {
            LOGGER.debug(String.format("grp(%d) has %d pure samples, perc(%.3f) of total(%d)",
                    mGroupId, mSampleCount, mSampleCount / (double) mAllSampleCount, mAllSampleCount));
        }
    }

    public void setRatioRangePercentile(double value ) { mRatioRangePercentile = value; }

    public SigOptimiser(int groupId, final List<SampleData> samples, final List<double[]> sampleCountsList, final double[] refRatios)
    {
        // assessment only rountine - how is this different??
        mGroupId = groupId;
        mIsValid = true;
        mHasChanged = false;
        mAssessmentOnly = true;

        if (sampleCountsList.isEmpty() || samples.size() != sampleCountsList.size())
        {
            mIsValid = false;
        }

        mBucketCount = refRatios.length;

        // no restriction to pure samples only
        mAllSampleCount = samples.size();

        mAllSamples = Lists.newArrayList();
        mAllSamples.addAll(samples);

        mSamples = Lists.newArrayList();
        mSamples.addAll(samples);
        mSampleCount = mSamples.size();

        mBucketInfo = new HashMap();

        mStartRatios = new double[mBucketCount];
        copyVector(refRatios, mStartRatios);
        mCurrentRatios = new double[mBucketCount];
        copyVector(refRatios, mCurrentRatios);
        mRangeMeanRatios = new double[mBucketCount];

        mNewBuckets = Lists.newArrayList();
        mRemovedBuckets = Lists.newArrayList();
        mInitialBuckets = Lists.newArrayList();
        mActiveBuckets = Lists.newArrayList();

        boolean addAllBuckets = sumVector(mStartRatios) == 0;

        for (int b = 0; b < mBucketCount; ++b)
        {
            if(addAllBuckets || mStartRatios[b] > 0)
                mInitialBuckets.add(b);
        }

        mActiveBuckets.addAll(mInitialBuckets);

        mCandidateNewBuckets = Lists.newArrayList();

        mLogVerbose = false;
        mCacheBucketInfo = true;

        // make a matrix from the sample counts
        mSampleCounts = new Matrix(mBucketCount, mAllSampleCount);
        mSampleNoiseRatio = new Matrix(mBucketCount, mAllSampleCount);

        mProposedAllocs = Lists.newArrayList();
        mProposedTotals = new double[mAllSampleCount];

        mZeroedSamples = Lists.newArrayList();

        mUnallocTotal = 0;
        mAllocTotal = 0;

        // capture the sample counts for all the current (non-zero ratio) buckets
        // the proposed totals are either zero if not yet determined, or whatever the bucket group has currently allocated per sample
        mSampleTotals = new double[mAllSampleCount];
        mUnallocBucketTotals = new double[mBucketCount];

        double[][] scData = mSampleCounts.getData();

        for (int s = 0; s < mAllSampleCount; ++s)
        {
            final double[] counts = sampleCountsList.get(s);

            for (Integer b : mActiveBuckets)
            {
                scData[b][s] = counts[b];
                mSampleTotals[s] += counts[b];
                mUnallocBucketTotals[b] += counts[b];
            }
        }

        mSampleCounts.cacheTranspose();

        mCountsTotal = sumVector(mSampleTotals);
        mRangesLow = new double[mBucketCount];
        mRangesHigh = new double[mBucketCount];
        mCalcRanges = new double[mBucketCount];
    }

    public final double[] getFittedRatios() { return mCurrentRatios; }
    public boolean isValid() { return mIsValid; }
    public boolean hasChanged() { return mHasChanged; }
    public List<Integer> getNewBuckets() { return mNewBuckets; }
    public List<Integer> getNonZeroBuckets() { return mActiveBuckets; }
    public double getAllocTotal() { return mAllocTotal; }
    public double[] getRevisedSampleAllocCounts(int sampleIndex) { return mProposedAllocs.get(sampleIndex); }
    public double getRevisedSampleAlloc(int sampleIndex) { return mProposedTotals[sampleIndex]; }
    public double[] getRatioRanges() { return mCalcRanges; }

    public void setLogVerbose(boolean toggle) { mLogVerbose = toggle; }
    public void setCacheBucketInfo(boolean toggle) { mCacheBucketInfo = toggle; }

    public final String getBucketInfo(int bucketId) { return mBucketInfo.get(bucketId); }

    public boolean optimiseBucketRatios(boolean calcOptimalRatios, boolean calcRangeDistributions)
    {
        if (!mIsValid)
            return false;

        LOGGER.debug(String.format("grp(%d) optimising for samples(%d of %d) buckets(%d) candidates(%d) count(%s alloc=%s unalloc=%s)",
                mGroupId, mSampleCount, mAllSampleCount, mActiveBuckets.size(), mCandidateNewBuckets.size(),
                sizeToStr(mCountsTotal), sizeToStr(mAllocTotal), sizeToStr(mUnallocTotal)));

        logStats(false);
        logRatios();

        if(mSampleCount < MIN_PURE_SAMPLE_COUNT)
            return true;

        removeLowCandidateBuckets();

        if(calcRangeDistributions)
        {
            calcBucketDistributions();
        }
        else
        {
            initVector(mRangesHigh, 0);
            initVector(mRangesLow, 0);
            initVector(mCalcRanges, 0);
            copyVector(mStartRatios, mRangeMeanRatios);
        }

        if(calcOptimalRatios)
            calcOptimalRatios();

        calcFinalRanges();

        if (!mIsValid)
            return false;

        if (mHasChanged)
        {
            updateSampleAllocations();

            logStats(false);
            logRatios();
        }

        return true;
    }

    private void updateSampleAllocations()
    {
        /*
        for (int s = 0; s < mAllSampleCount; ++s)
        {
            // fill in all missing details for this new bucket
            final SampleData sample = mAllSamples.get(s);
            double sampleCount = sample.getElevatedBucketCounts()[bucket];
            mSampleTotals[s] += sampleCount;
            scData[bucket][s] = sampleCount;

            // take the lower of the allocation or the actual
            double bucketRatio = mCurrentRatios[bucket];
            double optBucketAlloc = (rangeBoundMean * mProposedTotals[s]) / (1 - rangeBoundMean);
            double bucketAlloc = min(sampleCount, optBucketAlloc);
            mProposedTotals[s] += bucketAlloc;

            double[] sampleBucketAllocs = mProposedAllocs.get(s);
            sampleBucketAllocs[bucket] = bucketAlloc;
        }
        */
    }

    private void calcFinalRanges()
    {
        for (Integer b : mActiveBuckets)
        {
            if (mCalcRanges[b] < 0)
            {
                // to prevent negatives on the min ratio side
                mCalcRanges[b] = min(mCurrentRatios[b], 0.0005);

                // cap the range proportionally to the ratio
                mCalcRanges[b] = min(mCalcRanges[b], mCurrentRatios[b] * BUCKET_RANGE_MAX_PERCENT);
            }
        }
    }

    private void removeLowCandidateBuckets()
    {
        // remove any candidate bucket that is below a required min vs the total across all buckets to save on evaluating it
        double totalCount = sumVector(mUnallocBucketTotals);

        double maxRatio = 0;
        for(Integer bucket : mInitialBuckets)
        {
            maxRatio = max(maxRatio, mUnallocBucketTotals[bucket] / totalCount);
        }

        int index = 0;
        while(index < mCandidateNewBuckets.size())
        {
            int bucket = mCandidateNewBuckets.get(index);

            double bucketPerc = mUnallocBucketTotals[bucket] / totalCount;

            if(bucketPerc >= SMALL_RATIO_PERC_CUTOFF * maxRatio)
            {
                ++index;
                continue;
            }

            mCandidateNewBuckets.remove(index);
        }
    }

    public void optimiseRatiosAndRanges()
    {
        List<Integer> bucketIds = Lists.newArrayList();
        bucketIds.addAll(mActiveBuckets);
        bucketIds.addAll(mCandidateNewBuckets);

        mSampleSigContrib = new Matrix(mBucketCount, mSamples.size());
        double[][] sscData = mSampleSigContrib.getData();

        mSampleRatios = new Matrix(mBucketCount, mSamples.size());
        double[][] srData = mSampleRatios.getData();

        // for each sample, take the unalloc counts and convert them into ratios
        // then sort these ratios into ordered lists by bucket

        // first populate a new matrix with the unallocated counts for each sample in the relevant buckets
        for(int s = 0; s < mSamples.size(); ++s)
        {
            final SampleData sample = mSamples.get(s);
            final double[] counts = sample.getUnallocBucketCounts();

            double sampleBgTotal = 0; // total count across all buckets just in this sample

            for(Integer b : bucketIds)
            {
                double sbCount = counts[b];
                sscData[b][s] = sbCount;
                sampleBgTotal += sbCount;
            }

            for(Integer b : bucketIds)
            {
                double sbCount = counts[b];
                srData[b][s] = sbCount / sampleBgTotal;
            }
        }

        // now iteratively sort each bucket's ratios and cap any high values, then recalculate
        int iterations = 0;
        while(iterations < 10)
        {
            ++iterations;

            for (Integer b : bucketIds)
            {
                calcBucketRatioLimits(b);

                double ratioHighRange = mRangesHigh[b];

                // cap any sample contribution above the bucket high ratio and then recalc those sample's ratios
                int adjustedSampleCount = 0;
                for (int s = 0; s < mSamples.size(); ++s)
                {
                    double sampleRatio = srData[b][s];

                    if (sampleRatio > ratioHighRange)
                    {
                        ++adjustedSampleCount;
                        capSampleContributions(bucketIds, b, s, ratioHighRange);
                    }
                }
            }
        }
    }

    private void calcBucketRatioLimits(int bucket)
    {
        double[] bucketRatios = mSampleRatios.getCol(bucket);

        List<Integer> ratioSortedIndices = getSortedVectorIndices(bucketRatios, true);

        // examine the spread of ratios and take the lower Xth percentile as the cap
        int medianIndex = (int)floor(ratioSortedIndices.size() / 2);
        double rangeLowRatio = ratioSortedIndices.get(0);
        double medianRatio = ratioSortedIndices.get(medianIndex);
        double rangeHighRatio = 0;

        for(int i = medianIndex + 1; i < ratioSortedIndices.size(); ++i)
        {
            Integer samIndex = ratioSortedIndices.get(i);
            double ratio = bucketRatios[samIndex];

            if(ratio/medianRatio > DEFAULT_RATIO_RANGE_PERCENTILE)
                break;

            rangeHighRatio = ratio;
        }

        mRangesLow[bucket] = rangeLowRatio;
        mRangesHigh[bucket] = rangeHighRatio;
        mRangeMeanRatios[bucket] = medianRatio;

        if(mLogVerbose)
        {
            LOGGER.debug(String.format("grp(%d) bucket(%d %s) ratio(start=%.3f median=%.3f) span(%.3f -> %.3f)",
                    mGroupId, bucket, mStartRatios[bucket], medianRatio, rangeLowRatio, rangeHighRatio));
        }
    }

    private void capSampleContributions(List<Integer> bucketIds, int bucket, int samIndex, double ratioLimit)
    {
        double[][] sscData = mSampleSigContrib.getData();
        double[][] srData = mSampleRatios.getData();

        double sampleRatio = srData[bucket][samIndex];

        double sampleBgTotal = 0;

        for(Integer b1 : bucketIds)
        {
            sampleBgTotal += sscData[b1][samIndex];
        }

        double currentBucketCount = sscData[bucket][samIndex];
        double reducedBucketCount = floor((sampleBgTotal - currentBucketCount) / (1 - sampleRatio) * ratioLimit);
        sscData[bucket][samIndex] = reducedBucketCount;
        sampleBgTotal -= currentBucketCount - reducedBucketCount;

        // update sample's ratios using lowered bucket count
        for(Integer b1 : bucketIds)
        {
            double sbCount = sscData[b1][samIndex];
            srData[b1][samIndex] = sbCount / sampleBgTotal;
        }
    }

    private static int BUCKET_RATIO_SEGMENT_COUNT = 20; // for dividing ratio percents into small segments
    private static int BUCKET_RATIO_SEGMENT_OUTLIER_LOW = 0;
    private static int BUCKET_RATIO_SEGMENT_OUTLIER_HIGH = BUCKET_RATIO_SEGMENT_COUNT + 1;

    public void calcBucketDistributions()
    {
        int segmentCount = BUCKET_RATIO_SEGMENT_COUNT+2;
        mBucketRatioWeights = new Matrix(mBucketCount, segmentCount);
        mBucketRatioFrequencies = new Matrix(mBucketCount, segmentCount);
        Matrix sampleBucketRatios = new Matrix(mBucketCount, mSampleCount);
        Matrix sampleBucketWeights = new Matrix(mBucketCount, mSampleCount);

        List<Integer> bucketIds = Lists.newArrayList();
        bucketIds.addAll(mActiveBuckets);
        bucketIds.addAll(mCandidateNewBuckets);

        extractSampleRatioData(bucketIds, sampleBucketRatios, sampleBucketWeights);

        for (Integer b : bucketIds)
        {
            examineBucketDistribution(b, sampleBucketRatios.getRow(b), sampleBucketWeights.getRow(b));
        }

        // update and take the final ratios
        convertToPercentages(mCurrentRatios);
    }

    private void assignBucketRatioFrequencies(int bucket, double[] sampleRatios, double[] sampleWeights, double ratioLow, double ratioHigh)
    {
        // translate raw sample ratios into slots in a defined range

        double[][] brwData = mBucketRatioWeights.getData();
        double[][] brfData = mBucketRatioFrequencies.getData();

        double ratioSpan = ratioHigh - ratioLow;
        double segmentDistance = ratioSpan / BUCKET_RATIO_SEGMENT_COUNT;

        for(int s = 0; s < sampleRatios.length; ++s)
        {
            double sampleBucketRatio = sampleRatios[s];
            double sampleWeight = sampleWeights[s];
            int segmentIndex = 0;

            if(sampleBucketRatio < ratioLow)
            {
                segmentIndex = BUCKET_RATIO_SEGMENT_OUTLIER_LOW;
            }
            else if(sampleBucketRatio > ratioHigh)
            {
                segmentIndex = BUCKET_RATIO_SEGMENT_OUTLIER_HIGH;
            }
            else
            {
                // find the closest segment to records this info
                segmentIndex = (int)round((sampleBucketRatio - ratioLow) / segmentDistance) + 1;
            }

            ++brfData[bucket][segmentIndex];
            brwData[bucket][segmentIndex] += sampleWeight;
        }
    }

    private void extractSampleRatioData(final List<Integer> bucketIds, Matrix sampleBucketRatios, Matrix sampleBucketWeights)
    {
        double[][] ratioData = sampleBucketRatios.getData();
        double[][] weightData = sampleBucketWeights.getData();

        for(int samIndex = 0; samIndex < mSamples.size(); ++samIndex)
        {
            final SampleData sample = mSamples.get(samIndex);
            final double[] counts = mAssessmentOnly ? mSampleCounts.getCol(samIndex) : sample.getUnallocBucketCounts();

            // total count across all buckets just in this sample
            double sampleBgTotal = bucketIds.stream().mapToDouble(x -> counts[x]).sum();

            for(Integer b : bucketIds)
            {
                double sbCount = counts[b];
                double rawBucketRatio = sbCount / sampleBgTotal;
                double sampleWeight = sbCount / mUnallocBucketTotals[b];

                ratioData[b][samIndex] = rawBucketRatio;
                weightData[b][samIndex] = sampleWeight;
            }
        }
    }

    private void examineBucketDistribution(int bucket, double[] rawBucketRatios, double[] sampleBucketWeights)
    {
        // define percentile bounds (eg to capture the middle 60% of samples)
        double startPercentile = (1 - mRatioRangePercentile) * 0.5;
        double endPercentile = 1 - startPercentile;

        // repeat again but this time using the raw bucket ratios for accuracy
        List<Integer> ratioSortedIndices = getSortedVectorIndices(rawBucketRatios, true);

        double cumulativeWeight = 0;
        double rangeMinRatio = 0;
        double rangeMaxRatio = 0;
        double rangeWeightTotal = 0;
        double rangeRatioWeightTotal = 0;
        double medianRatio = 0;

        if(!doublesEqual(sumVector(sampleBucketWeights), 1.0))
        {
            LOGGER.warn("bg({}) has invalid sample bucket weights total({})", mGroupId, sumVector(sampleBucketWeights));
        }

        for(Integer samIndex : ratioSortedIndices)
        {
            double sampleBucketRatio = rawBucketRatios[samIndex];
            double sbWeight = sampleBucketWeights[samIndex];

            cumulativeWeight += sbWeight;

            if(cumulativeWeight < startPercentile)
                continue;

            if (rangeMinRatio == 0)
            {
                rangeMinRatio = sampleBucketRatio;
            }

            if(medianRatio == 0 && cumulativeWeight >= 0.5)
                medianRatio = sampleBucketRatio;

            if(cumulativeWeight < endPercentile)
            {
                rangeWeightTotal += sbWeight;
                rangeRatioWeightTotal += sbWeight;
            }
            else
            {
                if(rangeMaxRatio == 0)
                    rangeMaxRatio = sampleBucketRatio;
            }
        }

        double bucketTotal = mUnallocBucketTotals[bucket];

        if(rangeWeightTotal == 0)
        {
            LOGGER.warn(String.format("grp(%d) bucket(%d) couldn't calc mean: range(%.3f -> %.3f) total(%s)",
                    mGroupId, bucket, rangeMinRatio, rangeMaxRatio, sizeToStr(bucketTotal)));
            return;
        }

        double rangeBoundMean = rangeRatioWeightTotal / rangeWeightTotal;

        mRangesLow[bucket] = rangeMinRatio;
        mRangesHigh[bucket] = rangeMaxRatio;
        mRangeMeanRatios[bucket] = rangeBoundMean;

        double minRange = min(rangeBoundMean - rangeMinRatio, rangeMaxRatio - rangeBoundMean);

        boolean isCandidate = mCandidateNewBuckets.contains(bucket);

        double countsTotal = rangeWeightTotal * bucketTotal;

        if(mLogVerbose)
        {
            LOGGER.debug(String.format("grp(%d) bucket(%d %s) startRatio(%.3f) calc(avg=%.3f med=%.3f) span(%.3f -> %.3f act=%.3f perc=%.3f) total(%.3f, %s of %s)",
                    mGroupId, bucket, isCandidate ? "cand" : "init", mStartRatios[bucket], rangeBoundMean, medianRatio, rangeMinRatio, rangeMaxRatio, minRange, minRange/rangeBoundMean,
                    rangeWeightTotal, sizeToStr(countsTotal), sizeToStr(bucketTotal)));
        }

        if(mCacheBucketInfo)
        {
            assignBucketRatioFrequencies(bucket, rawBucketRatios, sampleBucketWeights, rangeMinRatio, rangeMaxRatio);

            double[][] brwData = mBucketRatioWeights.getData();
            double[][] brfData = mBucketRatioFrequencies.getData();

            // cache this data for external logging
            String bucketDataStr = "";

            // various stats
            bucketDataStr += String.format("%.0f,%.0f,%.4f,%.4f,%.4f,%.4f,%.4f",
                    bucketTotal, countsTotal, mStartRatios[bucket], rangeBoundMean, rangeMinRatio, rangeMaxRatio, minRange,
                    brwData[bucket][BUCKET_RATIO_SEGMENT_OUTLIER_LOW], brwData[bucket][BUCKET_RATIO_SEGMENT_OUTLIER_HIGH]);

            // log the actual frequencies for both counts and weights
            for(int i = 1; i <= BUCKET_RATIO_SEGMENT_COUNT; ++i)
            {
                bucketDataStr += String.format(",%.3f", brwData[bucket][i]);
            }

            for(int i = 1; i <= BUCKET_RATIO_SEGMENT_COUNT; ++i)
            {
                bucketDataStr += String.format(",%.0f", brfData[bucket][i]);
            }

            mBucketInfo.put(bucket, bucketDataStr);
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
        double rangePerc = (rangeMaxRatio - rangeMinRatio) / rangeBoundMean;
        double maxRatioPerc = (rangeBoundMean / mRefMaxRatio);

        boolean bucketValid = (maxRatioPerc >= SMALL_RATIO_PERC_CUTOFF) && (rangePerc <= BUCKET_RANGE_MAX_PERCENT * 2);

        if(!bucketValid)
            return;

        mNewBuckets.add(bucket);
        // mActiveBuckets.add(bucket); // for now candidates are not actually added
        mHasChanged = true;

        mCurrentRatios[bucket] = rangeBoundMean;
        mCalcRanges[bucket] = minRange;
    }

    public static final String getBucketInfoHeader()
    {
        String headers = "Total,Alloc,SigRatio,MeanRatio,MinRatio,MaxRatio,Range,LowOutliers,HighOutliers";

        for(int i = 0; i < BUCKET_RATIO_SEGMENT_COUNT; ++i)
        {
            headers += String.format(",SW_%d", i);
        }

        for(int i = 0; i < BUCKET_RATIO_SEGMENT_COUNT; ++i)
        {
            headers += String.format(",SF_%d", i);
        }

        return headers;
    }

    private void calcOptimalRatios()
    {
        // now that a mean and range has been established, test each bucket for different ratio values within that range
        double[] currentRatios = new double[mBucketCount];
        double[] bestRatios = new double[mBucketCount];
        copyVector(mCurrentRatios, currentRatios);
        copyVector(mCurrentRatios, bestRatios);
        double[] testSampleAllocs = new double[mAllSampleCount];

        int bucketIncr = 5;

        double bestAllocTotal = 0;
        List<Integer> bucketsToRemove = Lists.newArrayList();
        List<Double> testRatios = Lists.newArrayList();

        for (Integer bucket : mActiveBuckets)
        {
            double startRatio = bestRatios[bucket]; // start with the existing ratio

            double ratioMean = mRangeMeanRatios[bucket];
            double rangeMin = mRangesLow[bucket];
            double rangeMax = mRangesHigh[bucket];
            double range = max(ratioMean - rangeMin, ratioMean + rangeMax);

            boolean changed = false;
            double bestRatio = 0;
            double testMin = 1;
            double testMax = 0;
            double bucketBestAllocTotal = 0;
            double bucketBestRatio = 0;

            // test a range of values around the starting ratio
            // and repeat the process with a tighter range centered around the new best ratio the second time through
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

                    if (testRatio > 0 && testRatio >= rangeMin && testRatio <= rangeMax)
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
            if(mInitialBuckets.contains(bucket))
            {
                LOGGER.warn("grp({}) removing initial bucket({})", mGroupId, bucket);
                // continue;
            }

            mActiveBuckets.remove(bucket);
            mRemovedBuckets.add(bucket);

            for (int s = 0; s < mAllSampleCount; ++s)
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

    private double calcTotalAllocatedCounts(final double[] ratios, double[] allocatedSampleCounts, List<Integer> zeroedSamples)
    {
        // goes through each sample and applies the input ratios to calculate what can be allocated
        // return data:
        // - allocated bucket counts - per-bucket counts with the input ratios
        // - zeroed samples - any sample which cannot allocate anything with the input ratios
        // - total allocated count across all samples

        final double[][] scData = mSampleCounts.getData();
        initVector(allocatedSampleCounts, 0);

        // find the limiting bucket and alloc
        double totalAllocated = 0;

        for (int s = 0; s < mAllSampleCount; ++s)
        {
            double minAlloc = 0;

            for (Integer b : mActiveBuckets)
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
            for (Integer b : mActiveBuckets)
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
                mActiveBuckets.size(), mNewBuckets.size(), mRemovedBuckets.size(), mSampleCount, mZeroedSamples.size()));
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

    public static double[] optimiseSampleFit(
            final SampleData sample, int groupId, final List<Integer> bucketIds, final double[] bucketRatios,
            final double[] ratioRanges, final double[] startAllocCounts, boolean logVerbose)
    {
        // tweak the ratios within the specified ranges to find the optimal alloc for this sample
        double startAllocTotal = sumVector(startAllocCounts);

        int bucketCount = bucketRatios.length;

        final double[] unallocCounts = sample.getUnallocBucketCounts();
        double[] currentUnallocCounts = new double[bucketCount];
        copyVector(unallocCounts, currentUnallocCounts);

        // first back out the currently allocated counts from this sig
        for (Integer bucket : bucketIds)
        {
            currentUnallocCounts[bucket] = currentUnallocCounts[bucket] + startAllocCounts[bucket];
        }

        double[] currentRatios = new double[bucketCount];
        double[] bestRatios = new double[bucketCount];
        copyVector(bucketRatios, currentRatios);
        copyVector(bucketRatios, bestRatios);

        int bucketIncr = 5;

        double bestAllocTotal = startAllocTotal;
        List<Double> testRatios = Lists.newArrayList();

        for (Integer bucket : bucketIds)
        {
            double startRatio = bestRatios[bucket]; // start with the existing ratio

            double range = ratioRanges[bucket];

            boolean changed = false;
            double bestRatio = 0;
            double testMin = 1;
            double testMax = 0;
            double bucketBestAllocTotal = 0;
            double bucketBestRatio = 0;

            // test a range of values around the starting ratio
            // and repeat the process with a tighter range centered around the new best ratio the second time through
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

                double rangeFrag = min(range, startRatio) / bucketIncr;

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

                    // now test the effect of this ratio adjustment
                    double allocTotal = 0;
                    boolean initialised = false;
                    for (Integer b : bucketIds)
                    {
                        if (currentRatios[b] == 0)
                            continue;

                        double alloc = currentUnallocCounts[b] / currentRatios[b];

                        if (!initialised || alloc < allocTotal)
                        {
                            initialised = true;
                            allocTotal = alloc;
                        }
                    }

                    // double allocTotal = calcTotalAllocatedCounts(currentRatios, testSampleAllocs, null);

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

                if(logVerbose)
                {
                    LOGGER.debug(String.format("sample(%d) bucket(%d) best ratio(%.4f converted=%.4f) vs start(%.4f) testRange(%.4f -> %.4f)",
                            sample.Id, bucket, bestRatio, bestRatios[bucket], bucketRatios[bucket], testMin, testMax));
                }
            }
        }

        if (bestAllocTotal <= startAllocTotal)
            return null;

        double allocPercImprove = (bestAllocTotal - startAllocTotal) / sample.getElevatedCount();
        LOGGER.debug(String.format("sample(%d) grp(%d) ratio-optim increased alloc count(%s -> %s) perc(%.3f of %s)",
                sample.Id, groupId, doubleToStr(startAllocTotal), doubleToStr(bestAllocTotal), allocPercImprove, sizeToStr(sample.getElevatedCount())));

        // produce a final optimal set of alloc counts
        double[] additionalAllocCounts = new double[bucketCount];
        for (Integer b : bucketIds)
        {
            additionalAllocCounts[b] = bestAllocTotal * bestRatios[b];
            additionalAllocCounts[b] = max(additionalAllocCounts[b] - startAllocCounts[b], 0);
        }

        return additionalAllocCounts;
    }


    /*
    private void calcBucketRatioRanges()
    {
        double[][] nrData = mSampleNoiseRatio.getData();
        double[][] scData = mSampleCounts.getData();

        for (int s = 0; s < mAllSampleCount; ++s)
        {
            final SampleData sample = mAllSamples.get(s);

            if(!mSamples.contains(sample))
                continue;

            final double[] allCounts = sample.getElevatedBucketCounts();
            final double[] noise = sample.getNoiseCounts();

            for (Integer b : mActiveBuckets)
            {
                if (allCounts[b] > 0)
                    nrData[b][s] = noise[b] / allCounts[b];
            }
        }

        // calculate a percentage range for noise around each bucket
        final double[] bucketRatios = mCurrentRatios;

        for (int b = 0; b < mBucketCount; ++b)
        {
            if (!mActiveBuckets.contains(b))
                continue;

            // for now skip the candidate ones added
            if (mNewBuckets.contains(b))
                continue;

            double percentTotal = 0;
            double noisePercentTotal = 0;
            double countsPercentTotal = 0;

            for (int s = 0; s < mAllSampleCount; ++s)
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
                //LOGGER.debug(String.format("grp(%d) bucket(%d) ratio(%.4f) range(noise=%.4f counts=%.4f)",
                //        mGroupId, b, bucketRatios[b], avgNoiseRange, avgCountsRange));
            }
        }
    }

    */

    /*
    private void testCandidateExtraBuckets()
    {
        if (mCandidateNewBuckets.isEmpty())
            return;

        // some similar groups threw up additional buckets to consider including
        // walk through all samples and extract an implied ratio for each of these buckets
        double[][] scData = mSampleCounts.getData();

        boolean recalcRequired = false;
        double minRatioRangePerc = 0.2;

        for (Integer bucket : mCandidateNewBuckets)
        {
            double[] sampleCounts = new double[mAllSampleCount];
            double[] sampleRatios = new double[mAllSampleCount];

            double minSampleRatio = 0;
            double maxSampleRatio = 0;

            for (int s = 0; s < mAllSampleCount; ++s)
            {
                final SampleData sample = mAllSamples.get(s);

                if(!mSamples.contains(sample))
                    continue;

                double noise = max(sample.getNoiseCounts()[bucket] - sample.getAllocNoiseCounts()[bucket], 0);
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

            double percIncrease = (maxAllocTotal - mAllocTotal) / mAllocTotal;

            if (percIncrease >= MIN_IMPROVE_PERCENT)
            {
                LOGGER.debug(String.format("grp(%d) adding candidate bucket(%d) optimal ratio(%.4f) with alloc improve(%.3f %s prev=%s) ratios(min=%.4f max=%.4f lsq=%.4f med=%.4f avg=%.4f)",
                        mGroupId, bucket, selectedRatio, percIncrease, doubleToStr(maxAllocTotal), doubleToStr(mAllocTotal),
                        minSampleRatio, maxSampleRatio, lsqRatio, trialRatio, avgRatio));

                mHasChanged = true;
                recalcRequired = true;
                mActiveBuckets.add(bucket);
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
    */

}
