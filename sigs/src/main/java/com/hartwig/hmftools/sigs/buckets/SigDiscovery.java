package com.hartwig.hmftools.sigs.buckets;

import static java.lang.Math.pow;

import static com.hartwig.hmftools.common.utils.VectorUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;
import static com.hartwig.hmftools.sigs.buckets.BaConfig.DEFAULT_SIG_RATIO_RANGE_PERCENT;
import static com.hartwig.hmftools.sigs.buckets.BaConfig.MIN_DISCOVERY_SAMPLE_COUNT;
import static com.hartwig.hmftools.sigs.buckets.BaConfig.MIN_GROUP_ALLOC_PERCENT;
import static com.hartwig.hmftools.sigs.common.CommonUtils.getDiffList;
import static com.hartwig.hmftools.sigs.common.CommonUtils.getMatchingList;
import static com.hartwig.hmftools.common.sigs.DataUtils.sizeToStr;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Matrix;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SigDiscovery
{
    private final BucketAnalyser mAnalyser;
    private BaConfig mConfig;
    private List<SampleData> mSampleData;
    private Matrix mSampleCounts;
    private int mBucketCount;
    private int mSampleCount;
    private double mTotalCount;

    private static final Logger LOGGER = LogManager.getLogger(SigDiscovery.class);

    public SigDiscovery(BucketAnalyser analyser)
    {
        mAnalyser = analyser;
    }

    public void setInitialState(BaConfig config, List<SampleData> sampleData, Matrix sampleCounts)
    {
        mConfig = config;
        mSampleCounts = sampleCounts;
        mSampleData = sampleData;
        mBucketCount = mSampleCounts.Rows;
        mSampleCount = mSampleCounts.Cols;
        mTotalCount = mSampleCounts.sum();
    }

    public void scaleCountsByMutationalLoad(double[] sampleCounts)
    {
        if(mConfig.MutLoadWeightFactor == 1)
            return;
    }

    public void formBucketGroupsFromSamplePairs()
    {
        // create unique groups from any subset of 2 samples' buckets if they are a close match
        // samples aren't necessarily allocated, only groups are made for the time-being
        LOGGER.debug("forming bucket groups from sample-pair reduced buckets");

        String bgTag = "subset";

        double[] sc1 = new double[mBucketCount];
        double[] sc2 = new double[mBucketCount];

        double minSampleCount = MIN_DISCOVERY_SAMPLE_COUNT * mTotalCount;

        int groupsCreated = 0;

        for (int samIndex1 = 0; samIndex1 < mSampleCount; ++samIndex1)
        {
            SampleData sample1 = mSampleData.get(samIndex1);

            if(sample1.isExcluded())
                continue;

            if(!mAnalyser.getReassessSamples().isEmpty() && !mAnalyser.getReassessSamples().contains(samIndex1))
                continue;

            // mConfig.logSample(samIndex1);

            double reqSam1AllocPercent = BucketAnalyser.minAllocPercent(sample1, false);

            if(sample1.getUnallocPercent() < reqSam1AllocPercent)
                continue;

            if(sample1.getTotalCount() < minSampleCount)
                continue;

            final List<Integer> bl1 = sample1.getUnallocBuckets();

            if (bl1.isEmpty())
                continue;

            // record the top matching other sample - this will be used to create the top-allocating group
            double maxAllocaTotal = 0;
            int maxOtherSample = -1;
            double maxCss = 0;
            List<Double> maxCombinedCounts = Lists.newArrayList();
            List<Integer> maxSharedBuckets = Lists.newArrayList();

            for (int samIndex2 = samIndex1 + 1; samIndex2 < mSampleCount; ++samIndex2)
            {
                SampleData sample2 = mSampleData.get(samIndex2);

                if(sample2.isExcluded())
                    continue;

                // mConfig.logSample(samIndex1) && mConfig.logSample(samIndex2);

                double reqSam2AllocPercent = BucketAnalyser.minAllocPercent(sample2, false);

                if(sample2.getUnallocPercent() < reqSam2AllocPercent)
                    continue;

                if(sample2.getTotalCount() < minSampleCount)
                    continue;

                final List<Integer> bl2 = sample2.getUnallocBuckets();

                if (bl2.isEmpty())
                    continue;

                List<Integer> commonBuckets = getMatchingList(bl1, bl2);
                int commonBucketCount = commonBuckets.size();

                if (commonBucketCount < mConfig.MinBucketCountOverlap)
                    continue;

                sample1.populateBucketCountSubset(sc1, commonBuckets);
                sample2.populateBucketCountSubset(sc2, commonBuckets);

                final double[] sam1ElevCounts = sample1.getUnallocBucketCounts();
                final double[] sam2ElevCounts = sample2.getUnallocBucketCounts();

                double sam1ElevTotal = sumVector(sc1);
                double sam2ElevTotal = sumVector(sc2);
                double elevatedTotal = sam1ElevTotal + sam2ElevTotal;

                if(sam1ElevTotal/sample1.getElevatedCount() < reqSam1AllocPercent || sam2ElevTotal/sample2.getElevatedCount() < reqSam2AllocPercent)
                    continue;

                double bcCss = BucketAnalyser.calcSharedCSS(sc1, sc2);

                boolean addGroup = false;

                List<Integer> removedBuckets = Lists.newArrayList();

                if (bcCss >= mConfig.HighCssThreshold)
                {
                    addGroup = true;
                }
                else if (commonBucketCount > mConfig.MinBucketCountOverlap)
                {
                    // attempt to find a match using less overlapping buckets
                    double[] cssResults = new double[commonBuckets.size()];

                    for (int i = 0; i < commonBuckets.size(); ++i)
                    {
                        Integer testBucket = commonBuckets.get(i);

                        sc1[testBucket] = 0;
                        sc2[testBucket] = 0;

                        // run CSS on this reduced set of buckets
                        cssResults[i] = BucketAnalyser.calcSharedCSS(sc1, sc2);

                        // add the bucket back in ahead of the next test
                        sc1[testBucket] = sam1ElevCounts[testBucket];
                        sc2[testBucket] = sam2ElevCounts[testBucket];
                    }

                    List<Integer> sortedCssIndices = getSortedVectorIndices(cssResults, false);

                    // now remove each buckets one by one, with the ones having the largest effect to raise CSS first
                    for (int i = 0; i < sortedCssIndices.size(); ++i)
                    {
                        int commonBucketIndex = sortedCssIndices.get(i);
                        int testBucket = commonBuckets.get(commonBucketIndex);

                        sc1[testBucket] = 0;
                        sc2[testBucket] = 0;

                        sam1ElevTotal = sumVector(sc1);
                        sam2ElevTotal = sumVector(sc2);
                        elevatedTotal = sam1ElevTotal + sam2ElevTotal;

                        if(sam1ElevTotal/sample1.getElevatedCount() < reqSam1AllocPercent || sam2ElevTotal/sample2.getElevatedCount() < reqSam2AllocPercent)
                            break;

                        // if(elevatedTotal < allBucketsTotal * 0.5)
                        //    break;

                        // run CSS on this reduced set of buckets
                        bcCss = BucketAnalyser.calcSharedCSS(sc1, sc2);

                        removedBuckets.add(testBucket);

                        if (bcCss >= mConfig.HighCssThreshold)
                        {
                            addGroup = true;
                            break;
                        }

                        if(commonBucketCount - removedBuckets.size() <= mConfig.MinBucketCountOverlap)
                            break;
                    }
                }

                if(!addGroup)
                    continue;

                double newAllocTotal = elevatedTotal * pow(bcCss, 2);

                if(newAllocTotal <= maxAllocaTotal)
                    continue;

                maxAllocaTotal = newAllocTotal;
                maxOtherSample = samIndex2;
                maxCombinedCounts.clear();
                maxCss = bcCss;
                maxSharedBuckets.clear();

                for (Integer bucket : commonBuckets)
                {
                    if (removedBuckets.contains(bucket))
                        continue;

                    maxSharedBuckets.add(bucket);
                    maxCombinedCounts.add(sc1[bucket] + sc2[bucket]);
                }

                //                LOGGER.debug(String.format("sample(%d) new top match sample2(%d) with buckets(s1=%d and s2=%d common=%d matched=%d) counts(s1=%s s2=%s) css(%.4f)",
                //                        samIndex1, samIndex2, bl1.size(), bl2.size(), commonBuckets.size(), maxSharedBuckets.size(),
                //                        sizeToStr(sam1ElevTotal), sizeToStr(sam2ElevTotal), bcCss));
            }

            // now convert these matched buckets and their counts into bucket ratios
            if(maxAllocaTotal == 0)
                continue;

            // now create a group from the best allocation for this sample
            BucketGroup bucketGroup = new BucketGroup(mAnalyser.getNextBucketId());
            bucketGroup.setTag(bgTag);
            bucketGroup.addInitialSample(samIndex1);
            bucketGroup.addInitialSample(maxOtherSample);
            bucketGroup.addBuckets(maxSharedBuckets);

            double[] bucketRatios = new double[mBucketCount];

            double countsTotal = 0;
            for (Double count : maxCombinedCounts)
            {
                countsTotal += count;
            }

            for (int index = 0; index < maxSharedBuckets.size(); ++index)
            {
                Integer bucket = maxSharedBuckets.get(index);
                bucketRatios[bucket] = maxCombinedCounts.get(index) / countsTotal;
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
                LOGGER.debug(String.format("added bg(%d) samples(%d and %d) with buckets(%d) css(%.4f) allocCalcTotal(%s)",
                        bucketGroup.getId(), samIndex1, maxOtherSample, maxSharedBuckets.size(), maxCss, sizeToStr(maxAllocaTotal)));
            }

            mAnalyser.addBucketGroup(bucketGroup);
            ++groupsCreated;
        }

        if(groupsCreated == 0)
        {
            LOGGER.debug("no sample-pair subset bucket groups created");
            return;
        }

        LOGGER.debug("created {} sample-pair {} groups", groupsCreated, bgTag);
    }

    public void formExcessBucketGroups()
    {
        // logic: rather than look for similarity in counts, put together a prospective group of
        // samples if they have the same or similar elevated buckets
        // then refine/tweak these groups using the SigOptimiser
        LOGGER.debug("forming bucket groups from excess elevated buckets");

        String bgTag = "excess";

        double[] sc1 = new double[mBucketCount];
        double[] sc2 = new double[mBucketCount];

        double minSampleCount = MIN_DISCOVERY_SAMPLE_COUNT * mTotalCount;
        double minBucketOverlapPerc = 0.80;

        List<Integer> allSamples = Lists.newArrayList();

        List<BucketGroup> newGroups = Lists.newArrayList();

        for (int samIndex1 = 0; samIndex1 < mSampleCount; ++samIndex1)
        {
            SampleData sample1 = mSampleData.get(samIndex1);

            if(sample1.isExcluded())
                continue;

            // if(!mReassessSamples.isEmpty() && !mReassessSamples.contains(samIndex1))
            //     continue;

            // mConfig.logSample(samIndex1);

            if(sample1.getUnallocPercent() < MIN_GROUP_ALLOC_PERCENT)
                continue;

            if(sample1.getTotalCount() < minSampleCount)
                continue;

            final List<Integer> bl1 = sample1.getUnallocBuckets();

            if (bl1.isEmpty())
                continue;

            double reqSam1AllocPercent = BucketAnalyser.minAllocPercent(sample1, false);

            List<Integer> sam1CoveredBuckets = Lists.newArrayList();

            // for check groups just created
            boolean addedToGroup = false;
            for(BucketGroup bucketGroup : newGroups)
            {
                final List<Integer> commonBuckets = getMatchingList(bucketGroup.getBucketIds(), bl1);

                if(bucketGroup.hasSample(samIndex1))
                {
                    // added to a new group with another sample
                    addedToGroup = true;

                    for(Integer bucket : commonBuckets)
                    {
                        if(!sam1CoveredBuckets.contains(bucket))
                            sam1CoveredBuckets.add(bucket);
                    }

                    continue;
                }

                if(commonBuckets.size() < minBucketOverlapPerc * bucketGroup.getBucketCount())
                    continue;

                final List<Integer> alreadyCovered = getMatchingList(commonBuckets, sam1CoveredBuckets);

                if(alreadyCovered.size() > 0.75 * sam1CoveredBuckets.size())
                    continue;

                sample1.populateBucketCountSubset(sc1, commonBuckets);

                double sam1ElevTotal = sumVector(sc1);

                if(sam1ElevTotal/sample1.getElevatedCount() < reqSam1AllocPercent)
                    continue;

                bucketGroup.addSample(samIndex1, sc1);
                addedToGroup = true;

                if(!allSamples.contains(samIndex1))
                    allSamples.add(samIndex1);

                for(Integer bucket : commonBuckets)
                {
                    if(!sam1CoveredBuckets.contains(bucket))
                        sam1CoveredBuckets.add(bucket);
                }
            }

            // for now if a sample if already part of a group, don't look for further matches
            // if(addedToGroup)
            //     continue;

            for (int samIndex2 = samIndex1 + 1; samIndex2 < mSampleCount; ++samIndex2)
            {
                SampleData sample2 = mSampleData.get(samIndex2);

                if (sample2.isExcluded())
                    continue;

                // mConfig.logSample(samIndex1) && mConfig.logSample(samIndex2)

                if (sample2.getUnallocPercent() < MIN_GROUP_ALLOC_PERCENT)
                    continue;

                if (sample2.getTotalCount() < minSampleCount)
                    continue;

                final List<Integer> bl2 = sample2.getUnallocBuckets();

                if (bl2.isEmpty())
                    continue;

                double reqSam2AllocPercent = BucketAnalyser.minAllocPercent(sample2, false);

                List<Integer> commonBuckets = getMatchingList(bl1, bl2);
                int commonBucketCount = commonBuckets.size();

                if (commonBucketCount < mConfig.MinBucketCountOverlap)
                    continue;

                if (commonBucketCount < minBucketOverlapPerc * bl1.size() || commonBucketCount < minBucketOverlapPerc * bl2.size())
                    continue;

                sample1.populateBucketCountSubset(sc1, commonBuckets);
                sample2.populateBucketCountSubset(sc2, commonBuckets);

                double sam1ElevTotal = sumVector(sc1);
                double sam2ElevTotal = sumVector(sc2);

                if(sam1ElevTotal/sample1.getElevatedCount() < reqSam1AllocPercent || sam2ElevTotal/sample2.getElevatedCount() < reqSam2AllocPercent)
                    continue;

                final List<Integer> alreadyCovered = getMatchingList(commonBuckets, sam1CoveredBuckets);

                if(alreadyCovered.size() > 0.75 * sam1CoveredBuckets.size())
                    continue;

                // now create a group from the best allocation for this sample
                BucketGroup bucketGroup = new BucketGroup(mAnalyser.getNextBucketId());
                bucketGroup.setTag(bgTag);
                bucketGroup.addInitialSample(samIndex1);
                bucketGroup.addInitialSample(samIndex2);
                bucketGroup.addBuckets(commonBuckets);

                bucketGroup.addSample(samIndex1, sc1);
                bucketGroup.addSample(samIndex2, sc2);

                if(!allSamples.contains(samIndex1))
                    allSamples.add(samIndex1);

                if(!allSamples.contains(samIndex2))
                    allSamples.add(samIndex2);

                if (mConfig.LogVerbose)
                {
                    LOGGER.debug(String.format("added bg(%d) samples(%d and %d) with buckets(s1=%d s2=%d match=%d)",
                            bucketGroup.getId(), samIndex1, samIndex2, bl1.size(), bl2.size(), commonBucketCount));
                }

                for(Integer bucket : commonBuckets)
                {
                    if(!sam1CoveredBuckets.contains(bucket))
                        sam1CoveredBuckets.add(bucket);
                }

                newGroups.add(bucketGroup);
            }
        }

        if(newGroups.isEmpty())
        {
            LOGGER.debug("no excess unalloc bucket groups created");
            return;
        }

        LOGGER.debug("created {} overlap-unallocated bucket groups, samplesIncluded({})", newGroups.size(), allSamples.size());

        Map<Integer, Integer> candidateBucketCounts = new HashMap();
        List<Integer> candidateBuckets = Lists.newArrayList();
        List<SampleData> samplesList = Lists.newArrayList();

        int minSamplesCount = 5;
        int minBucketsCount = 5;

        for(BucketGroup bucketGroup : newGroups)
        {
            if(bucketGroup.getSampleCount() < minSamplesCount || bucketGroup.getBucketCount() < minBucketsCount)
                continue;

            samplesList.clear();
            candidateBucketCounts.clear();

            for(Integer sampleId : bucketGroup.getSampleIds())
            {
                final SampleData sample = mSampleData.get(sampleId);
                samplesList.add(sample);

                List<Integer> diffBuckets = getDiffList(sample.getElevatedBuckets(), bucketGroup.getBucketIds());

                for(Integer bucket : diffBuckets)
                {
                    Integer bucketCount = candidateBucketCounts.get(bucket);
                    if(bucketCount == null)
                        candidateBucketCounts.put(bucket, 1);
                    else
                        candidateBucketCounts.put(bucket, bucketCount + 1);
                }
            }

            double requiredBucketCount = bucketGroup.getSampleCount() * 0.5;
            candidateBuckets.clear();

            for(Map.Entry<Integer, Integer> entry : candidateBucketCounts.entrySet())
            {
                if(entry.getValue() >= requiredBucketCount)
                    candidateBuckets.add(entry.getKey());
            }

            SigOptimiser sigOptim = new SigOptimiser(bucketGroup.getId(), samplesList, null, bucketGroup.getBucketRatios(), candidateBuckets);
            sigOptim.setLogVerbose(false);

            boolean isValid = sigOptim.optimiseBucketRatios(true, false);

            if(isValid)
            {
                if(sigOptim.hasChanged())
                {
                    bucketGroup.setBucketRatios(sigOptim.getFittedRatios());

                    if(mConfig.UseRatioRanges)
                        bucketGroup.setRatioRangePerc(DEFAULT_SIG_RATIO_RANGE_PERCENT);

                    bucketGroup.getBucketIds().clear();
                    bucketGroup.getBucketIds().addAll(sigOptim.getNonZeroBuckets());
                }

                LOGGER.debug("bg({}) added with samples({}) buckets({}) and potential alloc({})",
                        bucketGroup.getId(), bucketGroup.getSampleCount(), bucketGroup.getBucketCount(), sizeToStr(sigOptim.getAllocTotal()));

                mAnalyser.addBucketGroup(bucketGroup);
            }
        }
    }





}
