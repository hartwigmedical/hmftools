package com.hartwig.hmftools.sig_analyser.buckets;

import static java.lang.Math.round;

import static com.hartwig.hmftools.sig_analyser.buckets.BaConfig.DEFAULT_SIG_RATIO_RANGE_PERCENT;
import static com.hartwig.hmftools.sig_analyser.buckets.BucketAnalyser.calcSharedCSS;
import static com.hartwig.hmftools.sig_analyser.buckets.BucketGroup.BG_TYPE_BACKGROUND;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.calcRangeValue;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.convertToPercentages;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.copyVector;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.doublesEqual;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.sizeToStr;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.sumVector;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.sig_analyser.common.SigMatrix;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class BackgroundSigDiscovery
{
    private final BaConfig mConfig;
    private final List<SampleData> mSampleData;
    private final HashMap<String, List<Integer>> mCancerSamplesMap;

    private Map<String,BucketGroup> mCancerBucketGroups;
    private Map<String,Double> mCancerMutLoadThresholds;
    private Map<Integer, Integer> mNoiseRangeMap;
    private SigMatrix mBackgroundCounts;

    private final int mBucketCount;
    private final List<Integer> mAllBuckets;
    private int mNextBucketId;


    private static final Logger LOGGER = LogManager.getLogger(BackgroundSigDiscovery.class);

    public BackgroundSigDiscovery(
            final BaConfig config, final List<SampleData> sampleData, final HashMap<String, List<Integer>> cancerSamplesMap,
            Map<Integer, Integer> noiseRangeMap)
    {
        mConfig = config;
        mSampleData = sampleData;
        mCancerSamplesMap = cancerSamplesMap;
        mNoiseRangeMap = noiseRangeMap;

        mCancerBucketGroups = Maps.newHashMap();
        mCancerMutLoadThresholds = Maps.newHashMap();
        mNextBucketId = 0;

        mBucketCount = sampleData.get(0).getBucketCounts().length;
        mBackgroundCounts = new SigMatrix(mBucketCount, sampleData.size());

        mAllBuckets = Lists.newArrayList();

        for (int i = 0; i < mBucketCount; ++i)
        {
            mAllBuckets.add(i);
        }
    }

    public SigMatrix getBackgroundCounts() { return mBackgroundCounts; }
    public Map<String,BucketGroup> getCancerBucketGroups() { return mCancerBucketGroups; }
    public List<BucketGroup> getBucketGroups() { return mCancerBucketGroups.values().stream().collect(Collectors.toList()); }

    public void createBackgroundSigs()
    {
        // for each cancer type:
        for (Map.Entry<String, List<Integer>> entry : mCancerSamplesMap.entrySet())
        {
            final String cancerType = entry.getKey();
            final List<Integer> sampleIds = entry.getValue();

            // adjust the mutational load threshold as required so a minimum number of samples are used
            calcMutationalLoad(cancerType, sampleIds);

            initialiseSampleCounts(cancerType, sampleIds);

            // find potential background sigs using a CSS comparison method
            List<BucketGroup> candidateGroups = findCandidateGroups(cancerType, sampleIds);

            // find the optimum background sig by testing total allocation across all samples
            BucketGroup topGroup = findTopAllocatedGroup(cancerType, sampleIds, candidateGroups);

            // allocate each sample to the background group
            allocateSamples(cancerType, sampleIds, topGroup);

            mCancerBucketGroups.put(cancerType, topGroup);
        }
    }

    private void initialiseSampleCounts(final String cancerType, final List<Integer> sampleIds)
    {
        // before calculating potential allocations for this background signature. determine possible noise for each bucket
        double mutLoadThreshold = mCancerMutLoadThresholds.get(cancerType);

        for (Integer sampleId : sampleIds)
        {
            SampleData sample = mSampleData.get(sampleId);
            double[] bucketCounts = sample.getBucketCounts();
            double[] noiseValues = new double[mBucketCount];
            double sampleTotal = sample.getTotalCount();

            for (int i = 0; i < mBucketCount; ++i)
            {
                int bucketValue;
                if(sampleTotal > mutLoadThreshold)
                {
                    bucketValue = (int)round(bucketCounts[i] / sampleTotal * mutLoadThreshold);
                }
                else
                {
                    bucketValue = (int)bucketCounts[i];
                }


                noiseValues[i] = calcRangeValue(mNoiseRangeMap, bucketValue);
            }

            sample.setElevatedBucketCounts(sample.getBucketCounts(), noiseValues);
        }
    }

    private double calcMutationalLoad(final String cancerType, final List<Integer> sampleIds)
    {
        // add the counts from every sample with mutational load below the threshold
        // and then infer a background signature (literally bucket ratios) from those only
        int samplesIncluded = 0;

        double mutLoadThreshold = mConfig.MutationalLoadCap;

        // check how many samples will be below the mutational load threshold
        double[] ctSampleTotals = new double[sampleIds.size()];
        int sampleIndex = 0;

        for (Integer sampleId : sampleIds)
        {
            SampleData sample = mSampleData.get(sampleId);

            double sampleTotal = sample.getTotalCount();
            ctSampleTotals[sampleIndex++] = sampleTotal;

            if (sampleTotal <= mutLoadThreshold)
                ++samplesIncluded;
        }

        if(samplesIncluded < 5)
        {
            // increase the ML threshold to take the lowest 5 samples
            List<Integer> sortedRatioIndices = getSortedVectorIndices(ctSampleTotals, true);
            mutLoadThreshold = ctSampleTotals[sortedRatioIndices.get(4)];

            LOGGER.info("cancerType({}) using higher ML threshold({}) with only {} samples",
                    cancerType, String.format("%.0f", mutLoadThreshold));
        }

        mCancerMutLoadThresholds.put(cancerType, mutLoadThreshold);
        return mutLoadThreshold;
    }

    private List<BucketGroup> findCandidateGroups(final String cancerType, final List<Integer> sampleIds)
    {
        // create unique groups from any subset of 2 samples' buckets if they are a close match
        // samples aren't necessarily allocated, only groups are made for the time-being
        final List<BucketGroup> candidateGroups = Lists.newArrayList();

        double mutLoadThreshold = mCancerMutLoadThresholds.get(cancerType);

        String bgTag = "background";

        for (int samIndex1 = 0; samIndex1 < sampleIds.size() - 1; ++samIndex1)
        {
            int sampleId1 = sampleIds.get(samIndex1);
            SampleData sample1 = mSampleData.get(sampleId1);

            if(sample1.getTotalCount() > mutLoadThreshold)
                continue;

            final double[] sc1 = sample1.getBucketCounts();

            // record the top matching other sample - this will be used to create the top-allocating group
            int maxOtherSample = -1;
            double maxCss = 0;

            for (int samIndex2 = samIndex1 + 1; samIndex2 < sampleIds.size(); ++samIndex2)
            {
                int sampleId2 = sampleIds.get(samIndex2);
                SampleData sample2 = mSampleData.get(sampleId2);

                if(sample2.getTotalCount() > mutLoadThreshold)
                    continue;

                final double[] sc2 = sample2.getBucketCounts();

                double bcCss = calcSharedCSS(sc1, sc2);

                if (bcCss < mConfig.HighCssThreshold || bcCss < maxCss)
                    continue;

                maxOtherSample = sampleId2;
                maxCss = bcCss;
            }

            if(maxCss == 0)
                continue;

            // now create a group from the best allocation for this sample
            BucketGroup bucketGroup = new BucketGroup(mNextBucketId++);
            bucketGroup.setTag(bgTag);
            bucketGroup.addInitialSample(samIndex1);
            bucketGroup.addInitialSample(maxOtherSample);
            bucketGroup.addBuckets(mAllBuckets);

            double[] bucketRatios = new double[mBucketCount];

            SampleData sample2 = mSampleData.get(maxOtherSample);
            double countsTotal = sample1.getTotalCount() + sample2.getTotalCount();

            for (int i = 0; i < mBucketCount; ++i)
            {
                double combinedCount = sc1[i] + sample2.getBucketCounts()[i];
                bucketRatios[i] = combinedCount / countsTotal;
            }

            if(!doublesEqual(sumVector(bucketRatios), 1))
            {
                convertToPercentages(bucketRatios);
            }

            bucketGroup.setBucketRatios(bucketRatios);

            if(!bucketGroup.isValid())
            {
                LOGGER.warn("bg({}) has invalid ratios");
                continue;
            }

            if(mConfig.UseRatioRanges)
                bucketGroup.setRatioRangePerc(DEFAULT_SIG_RATIO_RANGE_PERCENT);

            if(mConfig.LogVerbose)
            {
                LOGGER.debug(String.format("added bg(%d) samples(%d & %d) totals(%d & %d) css(%.4f)",
                        bucketGroup.getId(), samIndex1, maxOtherSample, sample1.getTotalCount(), sample2.getTotalCount(), maxCss));
            }

            candidateGroups.add(bucketGroup);
        }

        LOGGER.debug("cancerType({}) created {} candidate bucket groups", cancerType, candidateGroups.size());
        return candidateGroups;
    }

    private BucketGroup findTopAllocatedGroup(final String cancerType, final List<Integer> sampleIds, final List<BucketGroup> candidateGroups)
    {
        BucketGroup topGroup = null;

        double mutLoadThreshold = mCancerMutLoadThresholds.get(cancerType);

        // first clear all existing allocations of samples to groups and vice versa
        for (BucketGroup bucketGroup : candidateGroups)
        {
            double[] bgRatios = new double[mBucketCount];
            copyVector(bucketGroup.getBucketRatios(), bgRatios); // won't be recomputed as sample counts are added

            for (Integer sampleId : sampleIds)
            {
                final SampleData sample = mSampleData.get(sampleId);

                double[] allocCounts = new double[mBucketCount];
                double allocCountTotal = sample.getPotentialCounts(bgRatios, mAllBuckets, bucketGroup.getRatioRanges(), allocCounts);
                double allocPercent = allocCountTotal / sample.getTotalCount();

                // adjust allocation down below the ML threshold if required
                if(allocCountTotal > mutLoadThreshold)
                {
                    double reduceRatio = mutLoadThreshold / allocCountTotal;
                    allocCountTotal = mutLoadThreshold;
                    allocPercent = 1;

                    for (int i = 0; i < mBucketCount; ++i)
                    {
                        allocCounts[i] *= reduceRatio;
                    }
                }

                bucketGroup.addPotentialAllocation(allocCountTotal);
                bucketGroup.addPotentialAdjAllocation(allocCountTotal * allocPercent);

                // add sample counts so ratios can be adjusted after the top group is selected
                bucketGroup.addSample(sampleId, allocCounts);
            }

            if(topGroup == null || bucketGroup.getPotentialAdjAllocation() > topGroup.getPotentialAdjAllocation())
            {
                topGroup = bucketGroup;
            }
        }

        if(topGroup == null)
        {
            LOGGER.info("cancerType({}) no top group found", cancerType);
            return null;
        }

        double[] initialRatios = topGroup.getBucketRatios();

        topGroup.recalcBucketRatios(mConfig.MutLoadWeightFactor);

        double ratioChangeCss = calcSharedCSS(initialRatios, topGroup.getBucketRatios());

        LOGGER.info(String.format("cancerType(%s) top group with allocation(%.0f) ratioChangeCss(%.3f)",
                cancerType, topGroup.getPotentialAllocation(), ratioChangeCss));

        return topGroup;
    }

    private void allocateSamples(final String cancerType, final List<Integer> sampleIds, BucketGroup bucketGroup)
    {
        bucketGroup.setGroupType(BG_TYPE_BACKGROUND);
        bucketGroup.setCancerType(cancerType);

        double mutLoadThreshold = mCancerMutLoadThresholds.get(cancerType);

        bucketGroup.clearSamples();

        // first clear all existing allocations of samples to groups and vice versa
        final double[] bgRatios = bucketGroup.getBucketRatios();

        for (Integer sampleId : sampleIds)
        {
            final SampleData sample = mSampleData.get(sampleId);

            double[] allocCounts = new double[mBucketCount];
            double allocCountTotal = sample.getPotentialCounts(bgRatios, mAllBuckets, bucketGroup.getRatioRanges(), allocCounts);

            if(allocCountTotal > mutLoadThreshold)
            {
                double reduceRatio = mutLoadThreshold / allocCountTotal;

                for (int i = 0; i < mBucketCount; ++i)
                {
                    allocCounts[i] *= reduceRatio;
                }
            }

            allocCountTotal = sample.allocateBucketCounts(allocCounts);
            double allocPercent = allocCountTotal / sample.getTotalCount();

            bucketGroup.addSample(sample.Id, allocCounts);

            sample.addBucketGroup(bucketGroup, allocPercent);

            // cache allocated background counts
            mBackgroundCounts.setCol(sampleId, allocCounts);

            LOGGER.debug(String.format("sample(%d) added to background group(%d) alloc(%s perc=%.3f of %s)",
                    sampleId, bucketGroup.getId(), sizeToStr(allocCountTotal), allocPercent, sizeToStr(sample.getTotalCount())));
        }

        mBackgroundCounts.cacheTranspose();
    }

    /* OLD METHODS

        if(cmd.hasOption(BA_SAMPLE_CALC_DATA_FILE))
        {
            final String fileName = cmd.getOptionValue(BA_SAMPLE_CALC_DATA_FILE);
            GenericDataCollection bgSampleAllocations = GenericDataLoader.loadFile(fileName, GD_TYPE_STRING);
            List<List<String>> dataSet = bgSampleAllocations.getStringData();

            for(final List<String> sampleData : dataSet)
            {
                // first 2 fields are sampleId and sampleName, other fields are data related
                if(sampleData.size() < 3)
                {
                    mHasErrors = true;
                    break;
                }

                int sampleId = Integer.parseInt(sampleData.get(SCD_COL_SAMPLE_ID));
                double bgAllocation = Double.parseDouble(sampleData.get(SCD_COL_BG_ALLOC));
                mSampleBgAllocations.add(sampleId, bgAllocation);
            }

            LOGGER.debug("loaded {} sample calc data fields", dataSet.size());
            mLoadedSampleCalcData = true;
        }



    private void formBackgroundBucketGroups()
    {
        if(!mConfig.UseBackgroundCounts)
            return;

        for (Map.Entry<String, List<Integer>> entry : mCancerSamplesMap.entrySet())
        {
            final String cancerType = entry.getKey();
            List<Integer> sampleIds = entry.getValue();

            if(!mConfig.SpecificCancer.isEmpty() && !mConfig.SpecificCancer.equals(cancerType))
                continue;

            assignToBackgroundBucketGroups(cancerType, sampleIds);
        }
    }

    private void assignToBackgroundBucketGroups(final String cancerType, final List<Integer> sampleIds)
    {
        LOGGER.debug("cancerType({}) creating background groups for {} samples", cancerType, sampleIds.size());

        final double[] bgBucketRatios = mCancerTypeBucketRatiosMap.get(cancerType);

        if(bgBucketRatios == null)
            return;

        List<Integer> bucketIds = Lists.newArrayList();
        for (int i = 0; i < mBucketCount; ++i)
        {
            if(bgBucketRatios[i] > 0)
                bucketIds.add(i);
        }

        BucketGroup bucketGroup = new BucketGroup(mNextBucketId++);
        bucketGroup.addBuckets(bucketIds);
        bucketGroup.setCancerType(cancerType);
        bucketGroup.setGroupType(BG_TYPE_BACKGROUND);
        // bucketGroup.setTag(BG_TYPE_BACKGROUND);
        bucketGroup.setBucketRatios(bgBucketRatios);

        if(mConfig.UseRatioRanges)
        {
            bucketGroup.setRatioRangePerc(DEFAULT_SIG_RATIO_RANGE_PERCENT);
        }

        mBackgroundGroups.add(bucketGroup);

        // for now add all samples to a single group per cancer type
        for (int index = 0; index < sampleIds.size(); ++index)
        {
            int sampleId = sampleIds.get(index);
            SampleData sample = mSampleData.get(sampleId);

            if(sample.isExcluded())
                continue;

            double[] sampleCounts = mBackgroundCounts.getCol(sampleId);
            double bgAllocTotal = sumVector(sampleCounts);
            double allocPerc = bgAllocTotal / sample.getTotalCount();
            sample.addBucketGroup(bucketGroup, allocPerc);

            LOGGER.debug(String.format("sample(%d) added to background group(%d) alloc(%s perc=%.3f of %s)",
                    sampleId, bucketGroup.getId(), sizeToStr(bgAllocTotal), allocPerc,
                    sizeToStr(sample.getTotalCount())));

            bucketGroup.addSample(sampleId, sampleCounts);
        }

        LOGGER.debug(String.format("background group(%d) samples(%d) totalAlloc(%s)",
                bucketGroup.getId(), bucketGroup.getSampleCount(), sizeToStr(bucketGroup.getTotalCount())));
    }

    private void calcBucketMedianData()
    {
        mSampleTotals = new double[mSampleCount];

        for (int i = 0; i < mSampleCount; ++i)
        {
            mSampleTotals[i] = sumVector(mSampleCounts.getCol(i));
        }

        mTotalCount = sumVector(mSampleTotals);

        mCancerTypeBucketRatiosMap = Maps.newHashMap();

        double[][] scData = mSampleCounts.getData();

        for (Map.Entry<String, List<Integer>> entry : mCancerSamplesMap.entrySet())
        {
            final String cancerType = entry.getKey();
            final List<Integer> sampleIds = entry.getValue();

            if(sampleIds.size() < 5)
            {
                LOGGER.info("skipping cancerType({}) with only {} samples", cancerType, sampleIds.size());

                if(mConfig.UseBackgroundCounts)
                    mHasErrors = true;
            }

            // add the counts from every sample with mutational load below the threshold
            // and then infer a background signature (literally bucket ratios) from those only
            int samplesIncluded = 0;

            double mutLoadThreshold = mConfig.MutationalLoadCap;

            // check how many samples will be below the mutational load threshold
            double[] ctSampleTotals = new double[sampleIds.size()];

            for (int j = 0; j < sampleIds.size(); ++j)
            {
                int sampleId = sampleIds.get(j);

                double sampleTotal = mSampleTotals[sampleId];
                ctSampleTotals[j] = sampleTotal;

                if (sampleTotal <= mutLoadThreshold)
                    ++samplesIncluded;
            }

            if(samplesIncluded < 5)
            {
                // increase the ML threshold to take the lowest 5 samples
                List<Integer> sortedRatioIndices = getSortedVectorIndices(ctSampleTotals, true);
                mutLoadThreshold = ctSampleTotals[sortedRatioIndices.get(4)];

                LOGGER.info("cancerType({}) using higher ML threshold({}) with only {} samples", cancerType, String.format("%.0f", mutLoadThreshold));
            }

            samplesIncluded = 0;
            List<double[]> sampleBucketRatios = Lists.newArrayList();

            for (int j = 0; j < sampleIds.size(); ++j)
            {
                int sampleId = sampleIds.get(j);

                double sampleTotal = mSampleTotals[sampleId];

                if (sampleTotal > mutLoadThreshold)
                    continue;

                ++samplesIncluded;

                double[] bucketRatios = new double[mBucketCount];

                for (int i = 0; i < mBucketCount; ++i)
                {
                    bucketRatios[i] = scData[i][sampleId] / sampleTotal;
                }

                sampleBucketRatios.add(bucketRatios);
            }

            // now convert back to average counts
            LOGGER.debug("cancerType({}) has {} low mutational load samples", cancerType, samplesIncluded);

            int medianIndex = samplesIncluded / 2;

            double[] medianBucketRatios = new double[mBucketCount];

            for (int i = 0; i < mBucketCount; ++i)
            {
                double[] bucketRatios = new double[samplesIncluded];

                for (int j = 0; j < sampleBucketRatios.size(); ++j)
                {
                    double[] sampleRatios = sampleBucketRatios.get(j);
                    bucketRatios[j] = sampleRatios[i];
                }

                // sort these and then take the median
                List<Integer> sortedRatioIndices = getSortedVectorIndices(bucketRatios, true);
                int medianRatioIndex = sortedRatioIndices.get(medianIndex);
                double medianRatio = bucketRatios[medianRatioIndex];
                medianBucketRatios[i] = medianRatio;
            }

            // convert to a percent
            convertToPercentages(medianBucketRatios);
            mCancerTypeBucketRatiosMap.put(cancerType, medianBucketRatios);
        }
    }

    private void calcSampleBackgroundCounts()
    {
        if(!mSampleBgAllocations.isEmpty())
            return;

        LOGGER.debug("calculating sample background counts");

        for(int s = 0; s < mSampleCount; ++s)
        {
            mSampleBgAllocations.add(s, 0.0);
        }

        if(!mConfig.UseBackgroundCounts)
            return;

        final double[][] scData = mSampleCounts.getData();

        for (Map.Entry<String, double[]> entry : mCancerTypeBucketRatiosMap.entrySet())
        {
            final String type = entry.getKey();
            List<Integer> sampleIds = mCancerSamplesMap.get(type);

            final double[] medianRatios = entry.getValue();

            for (Integer sampleId : sampleIds)
            {
                double[] sampleBgCounts = new double[mBucketCount];

                for (int i = 0; i < mBucketCount; ++i)
                {
                    int expectedCount = (int) round(medianRatios[i] * min(mSampleTotals[sampleId], mConfig.MutationalLoadCap));
                    sampleBgCounts[i] = min(expectedCount, scData[i][sampleId]);
                }

                double optimalBgCount = calcBestFitWithinProbability(sampleId, medianRatios, sampleBgCounts, 0.99, 0.01);
                mSampleBgAllocations.set(sampleId, optimalBgCount);
            }
        }
    }



    */

}
