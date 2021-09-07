package com.hartwig.hmftools.sigs.buckets;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.sigs.NoiseCalcs.calcRangeValue;
import static com.hartwig.hmftools.common.sigs.SigUtils.convertToPercentages;
import static com.hartwig.hmftools.common.utils.VectorUtils.copyVector;
import static com.hartwig.hmftools.common.utils.VectorUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;
import static com.hartwig.hmftools.sigs.buckets.BaConfig.CANCER_TYPE_OTHER;
import static com.hartwig.hmftools.sigs.buckets.BaConfig.DEFAULT_SIG_RATIO_RANGE_PERCENT;
import static com.hartwig.hmftools.sigs.buckets.BaConfig.SIG_SIMILAR_CSS;
import static com.hartwig.hmftools.sigs.buckets.BucketGroup.BG_TYPE_BACKGROUND;
import static com.hartwig.hmftools.sigs.common.CommonUtils.SIG_LOGGER;
import static com.hartwig.hmftools.sigs.common.CssRoutines.CSSR_I1;
import static com.hartwig.hmftools.sigs.common.CssRoutines.CSSR_I2;
import static com.hartwig.hmftools.sigs.common.CssRoutines.CSSR_VAL;
import static com.hartwig.hmftools.sigs.common.CssRoutines.getTopCssPairs;
import static com.hartwig.hmftools.common.sigs.DataUtils.doublesEqual;
import static com.hartwig.hmftools.common.sigs.DataUtils.sizeToStr;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Matrix;

public class BackgroundSigDiscovery
{
    private final BaConfig mConfig;
    private final List<SampleData> mSampleData;
    private final Map<String,List<Integer>> mCancerSamplesMap;

    private Map<String,BucketGroup> mCancerBucketGroups;
    private Map<String,Double> mCancerMutLoadThresholds;
    private Map<Integer, Integer> mNoiseRangeMap;
    private Matrix mBackgroundCounts;

    private final int mBucketCount;
    private final List<Integer> mAllBuckets;
    private int mNextBucketId;

    private static final String BACKGROUND_SIG_TAG = "background";

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
        mBackgroundCounts = new Matrix(mBucketCount, sampleData.size());

        mAllBuckets = Lists.newArrayList();

        for (int i = 0; i < mBucketCount; ++i)
        {
            mAllBuckets.add(i);
        }
    }

    public Matrix getBackgroundCounts() { return mBackgroundCounts; }
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

            BucketGroup medianRatiosGroup = createMedianRatiosGroup(cancerType, sampleIds);
            candidateGroups.add(medianRatiosGroup);

            // find the optimum background sig by testing total allocation across all samples
            BucketGroup topGroup = findTopAllocatedGroup(cancerType, sampleIds, candidateGroups);

            if(topGroup == null)
                continue;

            // allocate each sample to the background group
            allocateSamples(cancerType, sampleIds, topGroup);

            mCancerBucketGroups.put(cancerType, topGroup);
        }

        if(mConfig.RunBackgroundAnalysis)
        {
            compareGroups();
            checkSampleVsOtherSigs();
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

            SIG_LOGGER.info("cancerType({}) using higher ML threshold({}) with only {} samples",
                    cancerType, String.format("%.0f", mutLoadThreshold), samplesIncluded);
        }

        mCancerMutLoadThresholds.put(cancerType, mutLoadThreshold);
        return mutLoadThreshold;
    }

    private BucketGroup createMedianRatiosGroup(final String cancerType, final List<Integer> sampleIds)
    {
        // take the median ratio in each bucket from all samples below the ML threshold
        double mutLoadThreshold = mCancerMutLoadThresholds.get(cancerType);

        List<double[]> sampleBucketRatios = Lists.newArrayList();

        for (Integer sampleId : sampleIds)
        {
            SampleData sample = mSampleData.get(sampleId);

            double sampleTotal = sample.getTotalCount();

            if(sampleTotal == 0)
                continue;

            double[] bucketCounts = sample.getBucketCounts();

            if (sampleTotal > mutLoadThreshold)
                continue;

            double[] bucketRatios = new double[mBucketCount];

            for (int i = 0; i < mBucketCount; ++i)
            {
                bucketRatios[i] = bucketCounts[i] / sampleTotal;
            }

            sampleBucketRatios.add(bucketRatios);
        }

        int medianIndex = sampleBucketRatios.size() / 2;

        double[] medianBucketRatios = new double[mBucketCount];

        for (int i = 0; i < mBucketCount; ++i)
        {
            double[] bucketRatios = new double[sampleBucketRatios.size()];

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

        // convert to a percent and then add as a group
        convertToPercentages(medianBucketRatios);

        BucketGroup bucketGroup = new BucketGroup(mNextBucketId++);
        bucketGroup.setTag("medianRatios");
        bucketGroup.addBuckets(mAllBuckets);
        bucketGroup.setBucketRatios(medianBucketRatios);

        if(mConfig.UseRatioRanges)
            bucketGroup.setRatioRangePerc(DEFAULT_SIG_RATIO_RANGE_PERCENT);

        return bucketGroup;
    }

    private List<BucketGroup> findCandidateGroups(final String cancerType, final List<Integer> sampleIds)
    {
        // create unique groups from any subset of 2 samples' buckets if they are a close match
        // samples aren't necessarily allocated, only groups are made for the time-being
        final List<BucketGroup> candidateGroups = Lists.newArrayList();

        double mutLoadThreshold = mCancerMutLoadThresholds.get(cancerType);

        double maxAnyCss = 0;

        for (int samIndex1 = 0; samIndex1 < sampleIds.size() - 1; ++samIndex1)
        {
            int sampleId1 = sampleIds.get(samIndex1);
            SampleData sample1 = mSampleData.get(sampleId1);

            if(sample1.getTotalCount() > mutLoadThreshold)
                continue;

            final double[] sc1 = sample1.getBucketCounts();

            // record the top matching other sample - this will be used to create the top-allocating group
            SampleData bestOtherSample = null;
            double maxSampleCss = 0;

            for (int samIndex2 = samIndex1 + 1; samIndex2 < sampleIds.size(); ++samIndex2)
            {
                int sampleId2 = sampleIds.get(samIndex2);
                SampleData sample2 = mSampleData.get(sampleId2);

                if(sample2.getTotalCount() > mutLoadThreshold)
                    continue;

                final double[] sc2 = sample2.getBucketCounts();

                double bcCss = BucketAnalyser.calcSharedCSS(sc1, sc2);
                maxAnyCss = max(maxAnyCss, bcCss);

                if (bcCss < mConfig.HighCssThreshold || bcCss < maxSampleCss)
                    continue;

                bestOtherSample = sample2;
                maxSampleCss = bcCss;
            }

            if(bestOtherSample == null)
                continue;

            // now create a group from the best allocation for this sample
            BucketGroup bucketGroup = new BucketGroup(mNextBucketId++);
            bucketGroup.setTag(String.format("samples_%d_%d", sample1.Id, bestOtherSample.Id));
            bucketGroup.addInitialSample(sample1.Id);
            bucketGroup.addInitialSample(bestOtherSample.Id);
            bucketGroup.addBuckets(mAllBuckets);

            double[] bucketRatios = new double[mBucketCount];

            double countsTotal = sample1.getTotalCount() + bestOtherSample.getTotalCount();

            for (int i = 0; i < mBucketCount; ++i)
            {
                double combinedBucketCount = sc1[i] + bestOtherSample.getBucketCounts()[i];
                bucketRatios[i] = combinedBucketCount / countsTotal;
            }

            if(!doublesEqual(sumVector(bucketRatios), 1))
            {
                convertToPercentages(bucketRatios);
            }

            bucketGroup.setBucketRatios(bucketRatios);

            if(!bucketGroup.isValid())
            {
                SIG_LOGGER.warn("bg({}) has invalid ratios");
                continue;
            }

            if(mConfig.UseRatioRanges)
                bucketGroup.setRatioRangePerc(DEFAULT_SIG_RATIO_RANGE_PERCENT);

            if(mConfig.LogVerbose)
            {
                SIG_LOGGER.debug(String.format("added bg(%d) samples(%d & %d) totals(%d & %d) css(%.4f)",
                        bucketGroup.getId(), sample1.Id, bestOtherSample.Id,
                        sample1.getTotalCount(), bestOtherSample.getTotalCount(), maxSampleCss));
            }

            candidateGroups.add(bucketGroup);
        }

        SIG_LOGGER.debug("cancerType({}) created {} candidate bucket groups, maxCss({})",
                cancerType, candidateGroups.size(), String.format("%.3f", maxAnyCss));

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

                boolean highMutLoadSample = allocCountTotal > mutLoadThreshold;

                // adjust allocation down below the ML threshold if required
                if(highMutLoadSample)
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

                if(!highMutLoadSample)
                {
                    // add sample counts so ratios can be adjusted after the top group is selected
                    bucketGroup.addSample(sampleId, allocCounts);
                }
            }

            if(topGroup == null || bucketGroup.getPotentialAdjAllocation() > topGroup.getPotentialAdjAllocation())
            {
                topGroup = bucketGroup;
            }
        }

        if(topGroup == null)
        {
            SIG_LOGGER.info("cancerType({}) no top group found", cancerType);
            return null;
        }

        double[] initialRatios = new double[mBucketCount];
        copyVector(topGroup.getBucketRatios(), initialRatios);

        topGroup.recalcBucketRatios(mConfig.MutLoadWeightFactor);

        double ratioChangeCss = BucketAnalyser.calcSharedCSS(initialRatios, topGroup.getBucketRatios());

        SIG_LOGGER.info(String.format("cancerType(%s) top group(%d:%s) with allocation(%.0f) ratioChangeCss(%.3f)",
                cancerType, topGroup.getId(), topGroup.getTag(), topGroup.getPotentialAllocation(), ratioChangeCss));

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

        double sampleCountTotal = 0;
        double groupAllocTotal = 0;

        double lmlSampleCountTotal = 0;
        double lmlGroupAllocTotal = 0;

        for (Integer sampleId : sampleIds)
        {
            final SampleData sample = mSampleData.get(sampleId);

            double[] allocCounts = new double[mBucketCount];
            double allocCountTotal = sample.getPotentialCounts(bgRatios, mAllBuckets, bucketGroup.getRatioRanges(), allocCounts);

            boolean highMutLoadSample = allocCountTotal > mutLoadThreshold;

            if(highMutLoadSample)
            {
                double reduceRatio = mutLoadThreshold / allocCountTotal;

                for (int i = 0; i < mBucketCount; ++i)
                {
                    allocCounts[i] *= reduceRatio;
                }
            }

            allocCountTotal = sample.allocateBucketCounts(allocCounts);
            double allocPercent = allocCountTotal / sample.getTotalCount();

            sampleCountTotal += min(sample.getTotalCount(), mutLoadThreshold);
            groupAllocTotal += allocCountTotal;

            if(sample.getTotalCount() <= mutLoadThreshold)
            {
                lmlSampleCountTotal += min(sample.getTotalCount(), mutLoadThreshold);
                lmlGroupAllocTotal += allocCountTotal;
            }

            bucketGroup.addSample(sample.Id, allocCounts);

            sample.addBucketGroup(bucketGroup, allocPercent);

            // cache allocated background counts
            mBackgroundCounts.setCol(sampleId, allocCounts);

            SIG_LOGGER.debug(String.format("sample(%d) %s added to background group(%d) alloc(%s perc=%.3f of %s)",
                    sampleId, highMutLoadSample ? "high" : "low",
                    bucketGroup.getId(), sizeToStr(allocCountTotal), allocPercent, sizeToStr(sample.getTotalCount())));
        }

        SIG_LOGGER.debug(String.format("cancerType(%s samples=%d) group(%d:%s) low-ML alloc(%s perc=%.3f of %s) all alloc(%s perc=%.3f of %s)",
                cancerType, sampleIds.size(), bucketGroup.getId(), bucketGroup.getTag(),
                sizeToStr(lmlGroupAllocTotal), lmlGroupAllocTotal/lmlSampleCountTotal, sizeToStr(lmlSampleCountTotal),
                sizeToStr(groupAllocTotal), groupAllocTotal/sampleCountTotal, sizeToStr(sampleCountTotal)));

        mBackgroundCounts.cacheTranspose();
    }

    private void compareGroups()
    {
        final List<BucketGroup> bucketGroups = getBucketGroups();

        if(bucketGroups.size() <= 1)
            return;

        Matrix sigMatrix = new Matrix(mBucketCount, bucketGroups.size());

        for(int i = 0; i < bucketGroups.size(); ++i)
        {
            sigMatrix.setCol(i, bucketGroups.get(i).getBucketRatios());
        }

        double internalSigCssThreshold = SIG_SIMILAR_CSS * 0.9;

        // first the internally generated ones
        List<double[]> cssResults = getTopCssPairs(
                sigMatrix, sigMatrix, internalSigCssThreshold, true, true, true, false);

        if (cssResults.isEmpty())
        {
            SIG_LOGGER.debug("no similar proposed sigs from bucket groups");
        }
        else
        {
            for (final double[] result : cssResults)
            {
                double css = result[CSSR_VAL];
                int sigId1 = (int)result[CSSR_I1];
                int sigId2 = (int)result[CSSR_I2];
                final BucketGroup bg1 = bucketGroups.get(sigId1);
                final BucketGroup bg2 = bucketGroups.get(sigId2);

                SIG_LOGGER.debug(String.format("background sig(%d:%s) matches sig(%d:%s) with css(%.6f)",
                        sigId1, bg1.getCancerType(), bg2.getId(), bg2.getCancerType(), css));
            }
        }
    }

    private void checkSampleVsOtherSigs()
    {
        double[] allocCounts = new double[mBucketCount];

        Map<String,Integer> misAllocPairingCounts = Maps.newHashMap();
        Map<String,Integer> misAllocCancerCounts = Maps.newHashMap();

        for(SampleData sample : mSampleData)
        {
            double mutLoadThreshold = mCancerMutLoadThresholds.containsKey(sample.cancerType()) ?
                    mCancerMutLoadThresholds.get(sample.cancerType()) : mCancerMutLoadThresholds.get(CANCER_TYPE_OTHER);

            if(sample.getTotalCount() > mutLoadThreshold)
                continue;

            double allocPerc = sample.getAllocPercent();
            double allocTotal = sample.getAllocatedCount();
            double allocThreshold = allocTotal * 1.1;

            BucketGroup maxOtherGroup = null;
            double maxOtherAlloc = 0;

            for(Map.Entry<String,BucketGroup> entry : mCancerBucketGroups.entrySet())
            {
                final String cancerType = entry.getKey();

                if(cancerType.equals(sample.cancerType()))
                    continue;

                final BucketGroup bucketGroup = entry.getValue();

                double sigAllocTotal = sample.getPotentialCounts(
                        bucketGroup.getBucketRatios(), mAllBuckets, bucketGroup.getRatioRanges(), allocCounts);

                if(sigAllocTotal > allocThreshold && sigAllocTotal > maxOtherAlloc)
                {
                    maxOtherAlloc = sigAllocTotal;
                    maxOtherGroup = bucketGroup;
                }
            }

            if(maxOtherGroup != null)
            {
                String misAllocTypes = sample.cancerType() + "_" + maxOtherGroup.getCancerType();

                Integer count = misAllocPairingCounts.get(misAllocTypes);
                if(count == null)
                    misAllocPairingCounts.put(misAllocTypes, 1);
                else
                    misAllocPairingCounts.put(misAllocTypes, count+1);

                count = misAllocCancerCounts.get(sample.cancerType());

                if(count == null)
                    misAllocCancerCounts.put(sample.cancerType(), 1);
                else
                    misAllocCancerCounts.put(sample.cancerType(), count+1);

                SIG_LOGGER.info(String.format("sample(%s) own cancerType(%s alloc=%.3f) less than otherCT(%s alloc=%.3f) of total(%s)",
                        sample.Id, sample.cancerType(), allocPerc,
                        maxOtherGroup.getCancerType(), maxOtherAlloc/sample.getTotalCount(), sizeToStr(sample.getTotalCount())));
            }
        }

        for(Map.Entry<String,Integer> entry : misAllocCancerCounts.entrySet())
        {
            final String cancerType = entry.getKey();

            int ctSampleCount = mCancerSamplesMap.containsKey(cancerType) ?
                    mCancerSamplesMap.get(cancerType).size() : mCancerSamplesMap.get(CANCER_TYPE_OTHER).size();

            int misAllocCount = entry.getValue();
            double misAllocPerc = misAllocCount/(double)ctSampleCount;

            if(misAllocPerc >= 0.05 && misAllocCount >= 5)
            {
                SIG_LOGGER.info(String.format("mis-allocation cancerType(%s) counts(%s) perc(%.3f of %d)",
                        cancerType, entry.getValue(), misAllocPerc, ctSampleCount));
            }
        }

        for(Map.Entry<String,Integer> entry : misAllocPairingCounts.entrySet())
        {
            final String misAllocTypes = entry.getKey();
            final String[] cancerTypes = misAllocTypes.split("_");
            final String initCancerType = cancerTypes[0];

            int ctSampleCount = mCancerSamplesMap.containsKey(initCancerType) ?
                    mCancerSamplesMap.get(initCancerType).size() : mCancerSamplesMap.get(CANCER_TYPE_OTHER).size();

            int misAllocCount = entry.getValue();
            double misAllocPerc = misAllocCount/(double)ctSampleCount;

            if(misAllocPerc >= 0.05 && misAllocCount >= 5)
            {
                SIG_LOGGER.info(String.format("mis-allocation combo(%s) counts(%s) perc(%.3f of %d)",
                        entry.getKey(), entry.getValue(), misAllocPerc, ctSampleCount));
            }
        }
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

            SIG_LOGGER.debug("loaded {} sample calc data fields", dataSet.size());
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
        SIG_LOGGER.debug("cancerType({}) creating background groups for {} samples", cancerType, sampleIds.size());

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

            SIG_LOGGER.debug(String.format("sample(%d) added to background group(%d) alloc(%s perc=%.3f of %s)",
                    sampleId, bucketGroup.getId(), sizeToStr(bgAllocTotal), allocPerc,
                    sizeToStr(sample.getTotalCount())));

            bucketGroup.addSample(sampleId, sampleCounts);
        }

        SIG_LOGGER.debug(String.format("background group(%d) samples(%d) totalAlloc(%s)",
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
                SIG_LOGGER.info("skipping cancerType({}) with only {} samples", cancerType, sampleIds.size());

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

                SIG_LOGGER.info("cancerType({}) using higher ML threshold({}) with only {} samples", cancerType, String.format("%.0f", mutLoadThreshold));
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
            SIG_LOGGER.debug("cancerType({}) has {} low mutational load samples", cancerType, samplesIncluded);

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

        SIG_LOGGER.debug("calculating sample background counts");

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
