package com.hartwig.hmftools.data_analyser.calcs;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.data_analyser.calcs.BucketAnalyser.MIN_GROUP_ALLOC_PERCENT_LOWER;
import static com.hartwig.hmftools.data_analyser.calcs.BucketAnalyser.SAMPLE_ALLOCATED_PERCENT;
import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.CSSR_I1;
import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.CSSR_I2;
import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.CSSR_VAL;
import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.getTopCssPairs;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.copyVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getMatchingList;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getNewFile;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sizeToStr;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sumVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sumVectors;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.vectorMultiply;
import static com.hartwig.hmftools.data_analyser.types.BucketGroup.BG_TYPE_MAJOR;
import static com.hartwig.hmftools.data_analyser.types.BucketGroup.BG_TYPE_MINOR;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.data_analyser.loaders.GenericDataLoader;
import com.hartwig.hmftools.data_analyser.types.BucketGroup;
import com.hartwig.hmftools.data_analyser.types.GenericDataCollection;
import com.hartwig.hmftools.data_analyser.types.NmfMatrix;
import com.hartwig.hmftools.data_analyser.types.SampleData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class BaReporter
{
    // shared state with Bucket Analyser
    private GenericDataCollection mDataCollection;
    private NmfMatrix mSampleCounts;
    private double[] mSampleTotals;
    private double mTotalCount;
    private int mBucketCount;
    private int mSampleCount;
    private NmfMatrix mBackgroundCounts;
    private NmfMatrix mElevatedCounts; // actual - expected, capped at zero
    private double mElevatedCount;
    private double mAllocatedCount;
    private double mBackgroundCount;
    private NmfMatrix mSampleBucketRatios;

    private NmfMatrix mProposedSigs;
    private List<Integer> mSigToBgMapping;
    private NmfMatrix mReferenceSigs;
    private boolean mUsingRefSigs;

    private List<SampleData> mSampleData;
    private GenericDataCollection mExtSampleData;
    private HashMap<String,Integer> mExtCategoriesMap;
    private HashMap<String, List<Integer>> mCancerSamplesMap;

    private List<BucketGroup> mFinalBucketGroups;
    private List<BucketGroup> mBackgroundGroups;

    BufferedWriter mBgRatiosFileWriter;

    // config
    private String mOutputDir;
    private String mOutputFileId;

    private static int COL_CANCER_TYPE = 1;
    private static int CATEGORY_COL_COUNT = 3;
    private static String EFFECT_TRUE = "TRUE";

    // internal state
    private static final Logger LOGGER = LogManager.getLogger(BaReporter.class);

    public BaReporter()
    {
        mUsingRefSigs = false;
    }

    public void setInitialState(
            GenericDataCollection dataCollection, final String outputDir, final String outputFileId,
            final NmfMatrix sampleCounts, final List<SampleData> sampleData,
            final GenericDataCollection extSampleData, final HashMap<String,Integer> extCategoriesMap,final HashMap<String, List<Integer>> cancerSamplesMap,
            final List<BucketGroup> finalBucketGroups, final List<BucketGroup> backgroundGroups)
    {
        mDataCollection = dataCollection;
        mOutputDir = outputDir;
        mOutputFileId = outputFileId;
        mSampleCounts = sampleCounts;

        mTotalCount = mSampleCounts.sum();
        mBucketCount = mSampleCounts.Rows;
        mSampleCount = mSampleCounts.Cols;

        mAllocatedCount = 0;

        mSampleData = sampleData;
        mExtSampleData = extSampleData;
        mExtCategoriesMap = extCategoriesMap;
        mCancerSamplesMap = cancerSamplesMap;

        mFinalBucketGroups = finalBucketGroups;
        mBackgroundGroups = backgroundGroups;
    }

    public void setPreRunState(double[] sampleTotals, final NmfMatrix backgroundCounts, final NmfMatrix elevatedCounts, double totalCount, double elevatedCount)
    {
        mSampleTotals = sampleTotals;
        mBackgroundCounts = backgroundCounts;
        mElevatedCounts = elevatedCounts;
        mTotalCount = totalCount;
        mElevatedCount = elevatedCount;
    }

    public void setFinalState(final NmfMatrix proposedSigs, final List<Integer> sigToBgMapping)
    {
        mProposedSigs = proposedSigs;
        mSigToBgMapping = sigToBgMapping;

        try
        {
            if (mBgRatiosFileWriter != null)
                mBgRatiosFileWriter.close();
        }
        catch(IOException e)
        {
            LOGGER.error("failed to close bucket groups file");
        }
    }

    public double getAllocatedCount() { return mAllocatedCount; }

    public final NmfMatrix getReferenceSigs() { return mReferenceSigs; }

    public void loadReferenceSigs(final String filename, boolean usingRefSigs)
    {
        GenericDataCollection dataCollection = GenericDataLoader.loadFile(filename);
        mReferenceSigs = DataUtils.createMatrixFromListData(dataCollection.getData());
        mReferenceSigs.cacheTranspose();
        mUsingRefSigs = usingRefSigs;
    }

    public void logOverallStats()
    {
        mAllocatedCount = 0;
        int fullyAllocated = 0;
        double noiseAllocated = 0;

        for(final SampleData sample : mSampleData)
        {
            mAllocatedCount += sample.getAllocatedCount();

            if(sample.getAllocPercent() >= SAMPLE_ALLOCATED_PERCENT)
                ++fullyAllocated;

            noiseAllocated += sample.getAllocNoise();
        }

        mBackgroundCount = 0;
        for(final BucketGroup bucketGroup : mBackgroundGroups)
        {
            mBackgroundCount += bucketGroup.getTotalCount();
        }

        double elevatedCount = mTotalCount - mBackgroundCount;

        LOGGER.debug(String.format("overall: samples(%d alloc=%d) groups(%d) counts: total(%s) background(%s perc=%.3f) elevated(%s perc=%.3f) alloc(%s perc=%.3f) noise(%s perc=%.3f)",
                mSampleCount, fullyAllocated, mFinalBucketGroups.size(), sizeToStr(mTotalCount), sizeToStr(mBackgroundCount), mBackgroundCount/mTotalCount,
                sizeToStr(elevatedCount), elevatedCount/mTotalCount, sizeToStr(mAllocatedCount), mAllocatedCount/mElevatedCount,
                sizeToStr(noiseAllocated), noiseAllocated/mElevatedCount));
    }

    public void postRunAnalysis()
    {
        tagMajorGroups();
        logBucketGroups(true);
        logSampleResults();
        logWorstAllocatedSamples();
        logSimilarSampleContribs();
        logSigReconstructions();
    }

    public void logSampleResults()
    {
        final List<String> categories = mExtSampleData.getFieldNames();

        int fullyAllocCount = 0;
        int partiallyAllocCount = 0;
        int noMatchCount = 0;
        double countAllocated = 0; // the actual number of variants

        for(int sampleId = 0; sampleId < mSampleCount; ++sampleId)
        {
            final SampleData sample = mSampleData.get(sampleId);

            if(sample.isExcluded())
                continue;

            final List<Integer> samBucketList = sample.getElevatedBuckets();

            final List<String> sampleData = sample.getCategoryData();
            if(sampleData == null || sampleData.isEmpty())
                continue;

            String effects = "";
            int effectsCount = 0;
            for(int c = CATEGORY_COL_COUNT; c < sampleData.size(); ++c)
            {
                if(!sampleData.get(c).isEmpty())
                {
                    ++effectsCount;

                    if (!effects.isEmpty())
                        effects += ",";

                    if(sampleData.get(c).equals(EFFECT_TRUE))
                        effects += categories.get(c);
                    else
                        effects += categories.get(c) + "=" + sampleData.get(c);
                }
            }

            final String cancerType = sampleData.get(COL_CANCER_TYPE);

            countAllocated += sample.getAllocPercent() * sample.getElevatedCount();

            double bgTotal = sumVector(mBackgroundCounts.getCol(sampleId));
            boolean fullyAllocated = (sample.getAllocPercent() >= SAMPLE_ALLOCATED_PERCENT);
            boolean partiallyAllocated = !fullyAllocated && (sample.getAllocPercent() > 0.1);

            if(fullyAllocated)
                ++fullyAllocCount;

            if(partiallyAllocated)
                ++partiallyAllocCount;

            if(fullyAllocated || partiallyAllocated)
            {
                double noiseAllocPerc = sample.getAllocNoise() / sample.getNoiseTotal();

                LOGGER.debug(String.format("sample(%d: %s) %s allocated: groups(%d) buckets(%d unalloc=%d) count(total=%s bg=%s elev=%s alloc=%.3f noise=%.2f of %s) cancer(%s) effects(%d: %s)",
                        sampleId, sample.getSampleName(), fullyAllocated ? "fully" : "partially",
                        sample.getElevBucketGroups().size(), samBucketList.size(), sample.getUnallocBuckets().size(),
                        sizeToStr(sample.getTotalCount()), sizeToStr(bgTotal),
                        sizeToStr(sample.getElevatedCount()), sample.getAllocPercent(), noiseAllocPerc, sizeToStr(sample.getNoiseTotal()),
                        cancerType, effectsCount, effects));
                continue;
            }

            LOGGER.debug("sample({}: {}) unallocated with no match: buckets({}) count({} bg={} total={}) cancer({}) effects({}: {})",
                    sampleId, sample.getSampleName(), samBucketList.size(), sizeToStr(sample.getElevatedCount()), sizeToStr(bgTotal),
                    sizeToStr(sample.getTotalCount()), cancerType, effectsCount, effects);

            ++noMatchCount;
        }

        double percAllocated = countAllocated / mElevatedCount;

        LOGGER.debug(String.format("sample summary: total(%d) alloc(%d, %.3f of %s) partial(%d) unalloc(%d)",
                mSampleCount, fullyAllocCount, percAllocated, sizeToStr(mElevatedCount),
                partiallyAllocCount, noMatchCount));

        logOverallStats();
    }

    public void logBucketGroups(boolean verbose)
    {
        List<BucketGroup> bgList = mFinalBucketGroups;

        // sort by score before logging - for now log in order of creation (by highest allocation)
        // Collections.sort(bgList);

        // log top groups
        int maxToLog = verbose ? bgList.size() : min(bgList.size(), 40);

        if(verbose)
        {
            LOGGER.debug("logging all bucket groups of total({})", bgList.size());
        }
        else
        {
            LOGGER.debug("logging top {} bucket groups of total({})", maxToLog, bgList.size());
        }

        double totalCount = mElevatedCount;

        for (int i = 0; i < maxToLog; ++i)
        {
            BucketGroup bucketGroup = bgList.get(i);

            double allocPercTotal = 0;
            for(int sampleId : bucketGroup.getSampleIds())
            {
                double allocPerc = bucketGroup.getSampleCount(sampleId) / mSampleTotals[sampleId];
                allocPercTotal += allocPerc;
            }

            double avgAllocPerc = allocPercTotal / bucketGroup.getSampleIds().size();

            if (verbose)
            {
                String bucketIdsStr = "";

                // log top 20 initial buckets, then any extra ones added through the broadending routine
                final List<Integer> initBuckets = bucketGroup.getInitialBucketIds();
                int added = 0;
                List<Integer> descBucketRatioIndices = getSortedVectorIndices(bucketGroup.getBucketRatios(), false);
                for(Integer bucket : descBucketRatioIndices)
                {
                    if(!initBuckets.contains(bucket))
                        continue;

                    if (!bucketIdsStr.isEmpty())
                        bucketIdsStr += ", ";

                    bucketIdsStr += bucket;

                    ++added;
                    if(added >= 20)
                        break;
                }

                if (!bucketGroup.getExtraBucketIds().isEmpty())
                {
                    bucketIdsStr += " extra=" + bucketGroup.getExtraBucketIds().toString();
                }

                double groupPerc = bucketGroup.getTotalCount() / totalCount;

                String linkData = "";
                if(!bucketGroup.getTag().isEmpty())
                {
                    linkData += String.format("tag=%s ", bucketGroup.getTag());
                }

                LOGGER.debug(String.format("rank %d: bg(%d) %scancer(%s) samples(%d) variants(avg=%s avgAllocPerc=%.3f total=%s perc=%.3f) buckets(%d: %s) effects(%s)",
                        i, bucketGroup.getId(), linkData, bucketGroup.getCancerType(), bucketGroup.getSampleIds().size(),
                        sizeToStr(bucketGroup.getAvgCount()), avgAllocPerc, sizeToStr(bucketGroup.getTotalCount()), groupPerc,
                        bucketGroup.getBucketIds().size(), bucketIdsStr, bucketGroup.getEffects()));
            }
            else
            {
                LOGGER.debug(String.format("rank %d: %s bg(%d) score(%.0f size=%d purity=%.2f) samples(%d) buckets(%d)",
                        i, bucketGroup.getId(), bucketGroup.calcScore(), bucketGroup.getSize(), bucketGroup.getPurity(),
                        bucketGroup.getSampleIds().size(), bucketGroup.getBucketIds().size()));
            }
        }
    }

    public void logWorstAllocatedSamples()
    {
        List<int[]> worstAllocatedSamples = Lists.newArrayList();
        double[] allUnallocatedByBucket = new double[mBucketCount];
        double[] sigUnallocatedByBucket = new double[mBucketCount]; // only for buckets more than 10% unallocated

        int SAMPLE_ID = 0;
        int UNALLOC_AMT = 1;

        int nonExcludedSampleCount = 0;

        for(final SampleData sample : mSampleData)
        {
            if(sample.isExcluded())
                continue;

            ++nonExcludedSampleCount;
            final double[] unallocCounts = sample.getUnallocBucketCounts();
            sumVectors(unallocCounts, allUnallocatedByBucket);

            for(int b = 0; b < mBucketCount; ++b)
            {
                if(unallocCounts[b] >= 0.1 * sample.getElevatedBucketCounts()[b])
                    sigUnallocatedByBucket[b] += unallocCounts[b];
            }

            if(sample.getAllocPercent() > 0.95)
                continue;

            int unallocTotal = (int)round(sample.getUnallocPercent() * sample.getElevatedCount());
            int worstIndex = 0;
            for(; worstIndex < worstAllocatedSamples.size(); ++worstIndex)
            {
                if(unallocTotal > worstAllocatedSamples.get(worstIndex)[UNALLOC_AMT])
                    break;
            }

            if(worstIndex <= 100)
            {
                int[] worstData = { sample.Id, unallocTotal };
                worstAllocatedSamples.add(worstIndex, worstData);
            }
        }

        double totalUnallocated = 0;

        // log worst 2% or top 20 samples
        int maxSamplesToLog = (int)round(max(nonExcludedSampleCount * 0.02, 20));
        int sampleCount = min(worstAllocatedSamples.size(), maxSamplesToLog);

        for(int worstIndex = 0; worstIndex < sampleCount; ++worstIndex)
        {
            final int[] worstData = worstAllocatedSamples.get(worstIndex);
            final SampleData sample = mSampleData.get(worstData[SAMPLE_ID]);
            double unallocTotal = worstData[UNALLOC_AMT];

            totalUnallocated += sample.getUnallocatedCount();

            LOGGER.debug(String.format("%d: worst sample(%d: %s) cancer(%s) unallocated(%s of %s, perc=%.3f) percOfTotal(%.4f) groupCount(%d)",
                    worstIndex, sample.Id, sample.getSampleName(), sample.getCancerType(), sizeToStr(unallocTotal),
                    sizeToStr(sample.getElevatedCount()), sample.getUnallocPercent(), unallocTotal/mElevatedCount, sample.getElevBucketGroups().size()));
        }

        double allUnallocTotal = mElevatedCount - mAllocatedCount;

        LOGGER.debug(String.format("worst %d (perc=%.3f) samples have unalloc(%s) ofAllUnalloc(perc=%.3f of %s) ofTotal(%.3f)",
                sampleCount, sampleCount/(double)nonExcludedSampleCount, sizeToStr(totalUnallocated), totalUnallocated/allUnallocTotal,
                sizeToStr(allUnallocTotal), totalUnallocated/mAllocatedCount));

        // report worst allocated buckets in a similar way
        List<Integer> worstBucketIndices = getSortedVectorIndices(sigUnallocatedByBucket, false);
        for(int worstIndex = 0; worstIndex < 10; ++worstIndex)
        {
            int worstBucket = worstBucketIndices.get(worstIndex);
            double allBucketTotal = allUnallocatedByBucket[worstBucket];
            double sigBucketTotal = sigUnallocatedByBucket[worstBucket];
            double totalBucketCount = sumVector(mElevatedCounts.getRow(worstBucket));

            LOGGER.debug(String.format("%d: worst bucket(%d) unallocated signif(%s of %s, perc=%.3f ofTotal=%.3f) bucket allUnlloc(%s)",
                    worstIndex, worstBucket, sizeToStr(sigBucketTotal), sizeToStr(totalBucketCount),
                    sigBucketTotal/totalBucketCount, sigBucketTotal/mElevatedCount, sizeToStr(allBucketTotal)));
        }
    }

    public void logSimilarSampleContribs()
    {
        // for for reoccurring patterns in sig allocations across samples
        // only look at contributions above X%
        int sigCount = mFinalBucketGroups.size();

        NmfMatrix contribMatrix = new NmfMatrix(sigCount, mSampleCount);
        double[][] contribData = contribMatrix.getData();

        double minSigContribPerc = MIN_GROUP_ALLOC_PERCENT_LOWER;

        int[] sampleSigCount = new int[mSampleCount];

        for(int bgIndex = 0; bgIndex < mFinalBucketGroups.size(); ++bgIndex)
        {
            final BucketGroup bucketGroup = mFinalBucketGroups.get(bgIndex);

            // exclude background groups from the comparison
            if(mBackgroundGroups.contains(bucketGroup))
                continue;

            final List<Integer> sampleIds = bucketGroup.getSampleIds();
            final List<Double> sampleContribs = bucketGroup.getSampleCountTotals();

            for(int samIndex = 0; samIndex < sampleIds.size(); ++samIndex)
            {
                Integer sampleId = sampleIds.get(samIndex);
                double sampleContrib = sampleContribs.get(samIndex);

                if(sampleContrib/mSampleTotals[sampleId] < minSigContribPerc)
                    continue;

                contribData[bgIndex][sampleId] = sampleContrib;
                ++sampleSigCount[sampleId];
            }
        }

        // exclude samples with too few samples
        for(int sample = 0; sample < mSampleCount; ++sample)
        {
            if(sampleSigCount[sample] < 2)
            {
                for(int sig = 0 ; sig < sigCount; ++sig)
                {
                    contribData[sig][sample] = 0;
                }
            }
        }

        contribMatrix.cacheTranspose();

        // now run CSS on every sample combination
        // first the internally generated ones
        double cssThreshold = 0.9999;

        List<double[]> cssResults = getTopCssPairs(
                contribMatrix, contribMatrix, cssThreshold, false, true, true, false);

        if (cssResults.isEmpty())
        {
            LOGGER.debug("no similar sample contribution ratios");
        }
        else
        {
            Map<String,Integer> groupComboMap = new HashMap();

            for (int i = 0; i < cssResults.size(); ++i)
            {
                final double[] result = cssResults.get(i);
                int samId1 = (int)result[CSSR_I1];
                int samId2 = (int)result[CSSR_I2];

                final SampleData sample1 = mSampleData.get(samId1);
                final SampleData sample2 = mSampleData.get(samId2);

                List<Integer> sam1BgList = sample1.getElevBucketGroups().stream().filter(x -> x.isBackground() == false).map(BucketGroup::getId).collect(Collectors
                        .toList());
                List<Integer> sam2BgList = sample2.getElevBucketGroups().stream().filter(x -> x.isBackground() == false).map(BucketGroup::getId).collect(Collectors.toList());

                final List<Integer> bgOverlaps = getMatchingList(sam1BgList, sam2BgList);

                Collections.sort(bgOverlaps);

                String bgIdsStr = "";
                for (Integer bgId : bgOverlaps)
                {
                    if (!bgIdsStr.isEmpty())
                        bgIdsStr += ";";

                    bgIdsStr += bgId;
                }

                Integer grpRepeatCount = groupComboMap.get(bgIdsStr);
                if(grpRepeatCount == null)
                    groupComboMap.put(bgIdsStr, 1);
                else
                    groupComboMap.put(bgIdsStr, grpRepeatCount + 1);

                /*
                if(i < 20)
                {
                    LOGGER.debug(String.format("samples(%d and %d) similar sig contribs css(%.4f) groups(s1=%d s2=%d match=%d) ids(%s)",
                            samId1, samId2, result[CSSR_VAL], sam1BgList.size(), sam2BgList.size(), bgOverlaps.size(), bgIdsStr));
                }
                */
            }

            for(Map.Entry<String,Integer> entry : groupComboMap.entrySet())
            {
                if(entry.getValue() < 5)
                    continue;

                LOGGER.debug("bgGroups({}) repeated {} times", entry.getKey(), entry.getValue());
            }
        }
    }

    private static double MAJOR_GROUP_ALLOC_PERC = 0.02;
    private static double MAJOR_GROUP_SAMPLE_PERC = 0.05;

    public void tagMajorGroups()
    {
        double reqSampleCount = mTotalCount * MAJOR_GROUP_ALLOC_PERC;
        double reqSamples = mSampleCount * MAJOR_GROUP_SAMPLE_PERC;

        for (final BucketGroup bucketGroup : mFinalBucketGroups)
        {
            if(bucketGroup.isBackground())
                continue;

            if(!bucketGroup.getGroupType().isEmpty())
                continue;

            if(bucketGroup.getTotalCount() >= reqSampleCount && bucketGroup.getSampleIds().size() >= reqSamples)
            {
                bucketGroup.setGroupType(BG_TYPE_MAJOR);
            }
            else
            {
                bucketGroup.setGroupType(BG_TYPE_MINOR);
            }
        }
    }

    public void logSigReconstructions()
    {
        if(mUsingRefSigs)
            return;

        LOGGER.debug("testing sig reconstruction");

        // test whether any sigs can be reconstructed from the others
        SigContribOptimiser sigOptim = new SigContribOptimiser(mBucketCount, false, 0.99);

        double[] testGroupRatios = new double[mBucketCount];
        double[] testGroupNoise = new double[mBucketCount]; // unused
        List<double[]> ratiosCollection = Lists.newArrayList();
        List<Integer> sigIds = Lists.newArrayList();
        List<Integer> testSigBuckets = Lists.newArrayList();
        double minRatioThreshold = 0.001;

        for (final BucketGroup testGroup : mFinalBucketGroups)
        {
            if (testGroup.isBackground())
                continue;

            copyVector(testGroup.getBucketRatios(), testGroupRatios);
            vectorMultiply(testGroupRatios, 10000); // to make the other group contributions be a percentage of total

            testSigBuckets.clear();
            ratiosCollection.clear();
            sigIds.clear();

            for (int b = 0; b < mBucketCount; ++b)
            {
                if (testGroupRatios[b] > 0)
                    testSigBuckets.add(b);
            }

            for (final BucketGroup otherGroup : mFinalBucketGroups)
            {
                if (otherGroup.isBackground() || otherGroup == testGroup)
                    continue;

                double[] otherGroupRatios = new double[mBucketCount];
                copyVector(otherGroup.getBucketRatios(), otherGroupRatios);

                boolean hasDiffBuckets = false;
                for (int b = 0; b < mBucketCount; ++b)
                {
                    if (otherGroupRatios[b] > 0 && !testSigBuckets.contains(b))
                    {
                        if (otherGroupRatios[b] < minRatioThreshold)
                        {
                            otherGroupRatios[b] = 0; // skip this bucket and continue on
                            continue;
                        }

                        hasDiffBuckets = true;
                        break;
                    }
                }

                if (hasDiffBuckets)
                    continue;

                ratiosCollection.add(otherGroupRatios);
                sigIds.add(otherGroup.getId());
            }

            if (ratiosCollection.size() < 2)
            {
                // LOGGER.debug(String.format("bg(%d) insufficient overlapping sigs for reconstruction", testGroup.getId()));
                continue;
            }

            sigOptim.initialise(testGroup.getId(), testGroupRatios, testGroupNoise, ratiosCollection, 0.01, 0);
            sigOptim.setSigIds(sigIds);
            sigOptim.setLogVerbose(false);

            boolean validCalc = sigOptim.fitToSample();

            if (!validCalc)
                continue;

            if (sigOptim.getAllocPerc() < 0.8)
            {
                LOGGER.debug(String.format("bg(%d) achieved low reconstruction from %d sigs to %.3f percent",
                        testGroup.getId(), sigOptim.contributingSigCount(), sigOptim.getAllocPerc()));
                continue;
            }

            LOGGER.debug(String.format("bg(%d) achieved reconstruction from %d sigs to %.3f percent:",
                    testGroup.getId(), sigOptim.contributingSigCount(), sigOptim.getAllocPerc()));

            final double[] sigContribs = sigOptim.getContribs();
            for (int sig = 0; sig < ratiosCollection.size(); ++sig)
            {
                if(sigContribs[sig] == 0)
                    continue;

                LOGGER.debug(String.format("bg(%d) from sigId(%d: %d) contrib(%.3f) percent", testGroup.getId(), sig, sigIds.get(sig), sigContribs[sig] / 100));
                testGroup.addGroupLinks(String.format("recon_%d_%.2f", sigIds.get(sig), sigContribs[sig] / 100));
            }
        }

        // additionally check for any sig whose buckets are wholly contained within another sig and has a large overlap in samples
        for (int bgIndex1 = 0; bgIndex1 < mFinalBucketGroups.size(); ++bgIndex1)
        {
            final BucketGroup testGroup = mFinalBucketGroups.get(bgIndex1);

            if (testGroup.isBackground())
                continue;

            final List<Integer> bucketIds = testGroup.getBucketIds();
            final List<Integer> sampleIds = testGroup.getSampleIds();

            for (int bgIndex2 = bgIndex1 + 1; bgIndex2 < mFinalBucketGroups.size(); ++bgIndex2)
            {
                final BucketGroup otherGroup = mFinalBucketGroups.get(bgIndex2);

                if(otherGroup.isBackground() || otherGroup == testGroup)
                    continue;

                final List<Integer> otherBucketIds = otherGroup.getBucketIds();

                final List<Integer> commonBuckets = getMatchingList(bucketIds, otherBucketIds);

                if(commonBuckets.size() < otherBucketIds.size())
                    continue;

                final List<Integer> otherSampleIds = otherGroup.getSampleIds();

                final List<Integer> commonSamples = getMatchingList(sampleIds, otherSampleIds);

                if(commonSamples.size() >= 0.25 * sampleIds.size() || commonSamples.size() >= 0.25 * otherSampleIds.size())
                {
                    LOGGER.debug(String.format("bg(%d) has possible correction bg(%d) buckets(bg1=%d bg2=%d) sample(bg1=%d bg2=%d match=%d)",
                            testGroup.getId(), otherGroup.getId(), bucketIds.size(), otherBucketIds.size(),
                            sampleIds.size(), otherSampleIds.size(), commonSamples.size()));

                    otherGroup.addGroupLinks(String.format("corr_bg_%d", testGroup.getId()));
                }
            }
        }
    }

    public void compareSignatures()
    {
        double sigCompareCss = 0.90;

        if(mProposedSigs != null && !mUsingRefSigs)
        {
            // first the internally generated ones
            List<double[]> cssResults = getTopCssPairs(mProposedSigs, mProposedSigs, sigCompareCss, true, true, true, false);

            if (cssResults.isEmpty())
            {
                LOGGER.debug("no similar proposed sigs from bucket groups");
            }
            else
            {
                for (final double[] result : cssResults)
                {
                    double css = result[CSSR_VAL];
                    int sigId1 = (int)result[CSSR_I1];
                    int sigId2 = (int)result[CSSR_I2];
                    final BucketGroup bg1 = mFinalBucketGroups.get(mSigToBgMapping.get(sigId1));
                    final BucketGroup bg2 = mFinalBucketGroups.get(mSigToBgMapping.get(sigId2));

                    // ignore reporting similarities in the BG groups
                    if(mBackgroundGroups.contains(bg1) && mBackgroundGroups.contains(bg2))
                        continue;

                    LOGGER.debug(String.format("proposed sig(%s bg=%d: ct=%s eff=%s samples=%d) matches sig(%d bg=%d: ct=%s eff=%s samples=%d) with css(%.4f)",
                            sigId1, bg1.getId(), bg1.getCancerType(), bg1.getEffects(), bg1.getSampleIds().size(),
                            sigId2, bg2.getId(), bg2.getCancerType(), bg2.getEffects(), bg2.getSampleIds().size(), css));

                    bg1.addGroupLinks(String.format("css_bg_%d_%.4f", bg2.getId(), css));
                    bg2.addGroupLinks(String.format("css_bg_%d_%.4f", bg1.getId(), css));
                }
            }
        }

        if(mReferenceSigs == null || mProposedSigs == null)
            return;

        List<double[]> cssResults = getTopCssPairs(mReferenceSigs, mProposedSigs, sigCompareCss, false, false);

        if (cssResults.isEmpty())
        {
            LOGGER.debug("no similar sigs between external ref and bucket groups");
        }
        else
        {
            for (final double[] result : cssResults)
            {
                double css = result[CSSR_VAL];
                int externalSigId = (int)result[CSSR_I1] + 1; // bumped up to correspond to convention of starting with 1
                int proposedSigId = (int)result[CSSR_I2];
                final BucketGroup bucketGroup = mFinalBucketGroups.get(mSigToBgMapping.get(proposedSigId));

                LOGGER.debug(String.format("external ref sig(%d) matches bg(%d: ct=%s eff=%s samples=%d purity=%.2f) with css(%.4f)",
                        externalSigId, bucketGroup.getId(), bucketGroup.getCancerType(), bucketGroup.getEffects(),
                        bucketGroup.getSampleIds().size(), bucketGroup.getPurity(), css));

                bucketGroup.addRefSig(String.format("sig_%d_css=%.4f", externalSigId, css));
            }
        }
    }

    public void assessSampleGroupAllocations(final BucketGroup bucketGroup)
    {
        final List<double[]> sampleGroupAllocCounts = bucketGroup.getSampleCounts();
        final List<Integer> sampleIds = bucketGroup.getSampleIds();
        final List<Integer> bucketIds = bucketGroup.getBucketIds();
        final double[] bucketRatios = bucketGroup.getBucketRatios();

        List<Double> fixedRatioSet = Lists.newArrayList(); // for now just a percent, range 0 - 1

        double ratioIncrement = 0.05;
        double startRatio = 0;
        for (int i = 0; i <= 20; ++i)
        {
            fixedRatioSet.add(startRatio + ratioIncrement * i);
        }

        // fixed ratios count x 3 values: fixed ratio, sample count with this ratio, unalloc total with this ratio
        double[][] unallocVsTotalRatioFreqSet = new double[fixedRatioSet.size()][3];
        int FIXED_RATIO = 0;
        int SAM_COUNT = 1;
        int UNALLOC = 2;

        // log all bucket unalloc counts
        try
        {
            if(mBgRatiosFileWriter == null)
            {
                mBgRatiosFileWriter = getNewFile(mOutputDir, mOutputFileId + "_ba_bkt_unalloc_ratios.csv");

                mBgRatiosFileWriter.write("BgId,SampleCount,Bucket,Ratio,TotalCount,AllocCount,UnallocCount");

                // all unallocated bucket counts will be written out
                for(Double ratio : fixedRatioSet)
                {
                    mBgRatiosFileWriter.write(String.format(",%.2f", ratio));
                }

                for(Double ratio : fixedRatioSet)
                {
                    mBgRatiosFileWriter.write(String.format(",%.2f", ratio));
                }

                mBgRatiosFileWriter.newLine();
            }

            BufferedWriter writer = mBgRatiosFileWriter;

            // for(Integer bucket : bucketIds)
            for(int bucket = 0; bucket < mBucketCount; ++bucket)
            {
                for(int i = 0; i < fixedRatioSet.size(); ++i)
                {
                    Double ratio = fixedRatioSet.get(i);
                    unallocVsTotalRatioFreqSet[i][FIXED_RATIO] = ratio;
                    unallocVsTotalRatioFreqSet[i][SAM_COUNT] = 0;
                    unallocVsTotalRatioFreqSet[i][UNALLOC] = 0;
                }

                double bucketTotal = 0;
                double groupBucketTotal = 0;
                double unallocBucketTotal = 0;

                for (int samIndex = 0; samIndex < sampleIds.size(); ++samIndex)
                {
                    final SampleData sample = mSampleData.get(sampleIds.get(samIndex));
                    double groupBucketCount = sampleGroupAllocCounts.get(samIndex)[bucket];
                    double sampleBucketCount = sample.getElevatedBucketCounts()[bucket];
                    double unallocCount = sample.getUnallocBucketCounts()[bucket];

                    bucketTotal += sampleBucketCount;
                    groupBucketTotal += groupBucketCount;

                    if(unallocCount == 0)
                    {
                        unallocVsTotalRatioFreqSet[0][SAM_COUNT] += 1;
                        continue;
                    }

                    unallocBucketTotal += unallocCount;

                    // form a ratio from what has been allocated to this group vs what is left over
                    double unallocVsTotalRatio = unallocCount / sampleBucketCount;

                    // insert these into the frequency set
                    int index = 0;
                    while (index < fixedRatioSet.size() - 1)
                    {
                        if (abs(unallocVsTotalRatioFreqSet[index][0] - unallocVsTotalRatio) <= 0.5 * ratioIncrement)
                            break;

                        ++index;
                    }

                    unallocVsTotalRatioFreqSet[index][SAM_COUNT] += 1;
                    unallocVsTotalRatioFreqSet[index][UNALLOC] += unallocCount;
                }

                writer.write(String.format("%d,%d,%d,%.4f,%.0f,%.0f,%.0f",
                        bucketGroup.getId(), sampleIds.size(), bucket, bucketRatios[bucket], bucketTotal, groupBucketTotal, unallocBucketTotal));

                for(int i = 0; i < fixedRatioSet.size(); ++i)
                {
                    writer.write(String.format(",%.0f", unallocVsTotalRatioFreqSet[i][SAM_COUNT]));
                }

                for(int i = 0; i < fixedRatioSet.size(); ++i)
                {
                    writer.write(String.format(",%.0f", unallocVsTotalRatioFreqSet[i][UNALLOC]));
                }

                writer.newLine();
            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to unalloc bucket-ratios file: {}", e.toString());
        }
    }


    public void calcBucketDistributions()
    {
        double percIncrement = 0.02;
        int ratioSets = (int)round(1 / percIncrement);

        mSampleBucketRatios = new NmfMatrix(mBucketCount, ratioSets);
        double[][] sbrData = mSampleBucketRatios.getData();

        double[] bucketTotals = new double[mBucketCount];

        for(int b = 0; b < mBucketCount; ++b)
        {
            bucketTotals[b] = sumVector(mSampleCounts.getRow(b));
        }

        for(final SampleData sample : mSampleData)
        {
            double sampleTotal = sample.getTotalCount();

            final double[] counts = sample.getUnallocBucketCounts();
            // final double[] sampleCounts = sample.getElevatedBucketCounts();
            final double[] noise = sample.getCountRanges();
            final double[] allocNoise = sample.getAllocNoiseCounts();
            final List<Integer> elevBuckets = sample.getElevatedBuckets();

            for(int b = 0; b < mBucketCount; ++b)
            {
                if(!elevBuckets.contains(b))
                    continue;

                double sbCount = counts[b];

                if(sbCount == 0)
                    sbCount = max(noise[b] - allocNoise[b], 0);

                double rawBucketRatio = sbCount / sampleTotal;
                int ratioIndex = (int)round(rawBucketRatio / percIncrement);

                if(ratioIndex < 0 || ratioIndex > ratioSets)
                {
                    LOGGER.error(String.format("sample(%d) bucket(%d) invalid ratioIndex(%d) from sbCount(%.1f) and sampleTotal(%.1f)",
                            sample.Id, b, ratioIndex, sbCount, sampleTotal));
                    continue;
                }

                double sampleWeight = sbCount / bucketTotals[b];

                // double bucketRatio = round(rawBucketRatio / percIncrement) * percIncrement;

                sbrData[b][ratioIndex] += sampleWeight;
            }
        }
    }

}
