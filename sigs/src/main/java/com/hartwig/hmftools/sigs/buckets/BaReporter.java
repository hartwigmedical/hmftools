package com.hartwig.hmftools.sigs.buckets;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.stats.CosineSimilarity.calcCosineSim;
import static com.hartwig.hmftools.common.utils.MatrixUtils.createMatrixFromListData;
import static com.hartwig.hmftools.common.utils.VectorUtils.copyVector;
import static com.hartwig.hmftools.common.utils.VectorUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVectors;
import static com.hartwig.hmftools.common.utils.VectorUtils.vectorMultiply;
import static com.hartwig.hmftools.sigs.buckets.BaConfig.MIN_GROUP_ALLOC_PERCENT;
import static com.hartwig.hmftools.sigs.buckets.BaConfig.MIN_GROUP_ALLOC_PERCENT_LOWER;
import static com.hartwig.hmftools.sigs.buckets.BaConfig.SAMPLE_ALLOCATED_PERCENT;
import static com.hartwig.hmftools.sigs.buckets.BaConfig.SIG_SIMILAR_CSS;
import static com.hartwig.hmftools.sigs.common.CssRoutines.CSSR_I1;
import static com.hartwig.hmftools.sigs.common.CssRoutines.CSSR_I2;
import static com.hartwig.hmftools.sigs.common.CssRoutines.CSSR_VAL;
import static com.hartwig.hmftools.sigs.common.CssRoutines.getTopCssPairs;
import static com.hartwig.hmftools.sigs.common.CommonUtils.getMatchingList;
import static com.hartwig.hmftools.common.sigs.DataUtils.sizeToStr;
import static com.hartwig.hmftools.sigs.buckets.BucketGroup.BG_TYPE_MAJOR;
import static com.hartwig.hmftools.sigs.buckets.BucketGroup.BG_TYPE_MINOR;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.GenericDataLoader;
import com.hartwig.hmftools.common.utils.GenericDataCollection;
import com.hartwig.hmftools.common.utils.Matrix;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class BaReporter
{
    final BucketAnalyser mAnalyser;

    // shared state with Bucket Analyser
    private GenericDataCollection mDataCollection;
    private Matrix mSampleCounts;
    private double[] mSampleTotals;
    private double mTotalCount;
    private int mBucketCount;
    private int mSampleCount;
    private int mActiveSampleCount;
    private Matrix mBackgroundCounts;
    private Matrix mElevatedCounts; // actual - expected, capped at zero
    private double mElevatedCount;

    private double mTotalAllocatedCount;
    private Matrix mSampleBucketRatios;

    private Matrix mProposedSigs;
    private List<Integer> mSigToBgMapping;
    private Matrix mReferenceSigs;
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

    public BaReporter(final BucketAnalyser analyser)
    {
        mAnalyser = analyser;
        mUsingRefSigs = false;
    }

    public void setInitialState(
            GenericDataCollection dataCollection, final String outputDir, final String outputFileId,
            final Matrix sampleCounts, final List<SampleData> sampleData,
            final GenericDataCollection extSampleData, final HashMap<String,Integer> extCategoriesMap,
            final HashMap<String, List<Integer>> cancerSamplesMap,
            final List<BucketGroup> finalBucketGroups, final List<BucketGroup> backgroundGroups)
    {
        mDataCollection = dataCollection;
        mOutputDir = outputDir;
        mOutputFileId = outputFileId;
        mSampleCounts = sampleCounts;

        mTotalCount = mSampleCounts.sum();
        mBucketCount = mSampleCounts.Rows;
        mSampleCount = mSampleCounts.Cols;
        mActiveSampleCount = mSampleCount;

        mTotalAllocatedCount = 0;

        mSampleData = sampleData;
        mExtSampleData = extSampleData;
        mExtCategoriesMap = extCategoriesMap;
        mCancerSamplesMap = cancerSamplesMap;

        mFinalBucketGroups = finalBucketGroups;
        mBackgroundGroups = backgroundGroups;
    }

    public void setPreRunState(double[] sampleTotals, final Matrix backgroundCounts, final Matrix elevatedCounts, double totalCount, double elevatedCount, int activeSampleCount)
    {
        mSampleTotals = sampleTotals;
        mBackgroundCounts = backgroundCounts;
        mElevatedCounts = elevatedCounts;
        mTotalCount = totalCount;
        mElevatedCount = elevatedCount;
        mActiveSampleCount = activeSampleCount;
    }

    public void setFinalState(final Matrix proposedSigs, final List<Integer> sigToBgMapping)
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

    public double getTotalAllocatedCount() { return mTotalAllocatedCount; }

    public final Matrix getReferenceSigs() { return mReferenceSigs; }

    public void loadReferenceSigs(final String filename, boolean usingRefSigs)
    {
        GenericDataCollection dataCollection = GenericDataLoader.loadFile(filename);
        mReferenceSigs = createMatrixFromListData(dataCollection.getData());
        mReferenceSigs.cacheTranspose();
        mUsingRefSigs = usingRefSigs;
    }

    public void logOverallStats()
    {
        mTotalAllocatedCount = 0;
        int fullyAllocated = 0;
        double noiseAllocated = 0;
        double elevatedAllocatedCount = 0;
        int samplesUsingElevatedForAllocated = 0;

        for(final SampleData sample : mSampleData)
        {
            if(sample.usingElevatedForAllocation())
            {
                ++samplesUsingElevatedForAllocated;
            }

            elevatedAllocatedCount += sample.getAllocatedCount();

            if(sample.getAllocPercent() >= SAMPLE_ALLOCATED_PERCENT)
                ++fullyAllocated;

            noiseAllocated += sample.getAllocNoise();
        }

        double backgroundAlloc = 0;
        for(final BucketGroup bucketGroup : mBackgroundGroups)
        {
            backgroundAlloc += bucketGroup.getTotalCount();
        }

        if(samplesUsingElevatedForAllocated == mSampleData.size())
        {
            mTotalAllocatedCount = elevatedAllocatedCount + backgroundAlloc;
        }
        else
        {
            mTotalAllocatedCount = elevatedAllocatedCount;
        }

        double initBackgroundCount = mTotalCount - mElevatedCount;

        LOGGER.debug(String.format("overall: samples(%d alloc=%d) groups(%d) counts: total(%s) initial(bg=%s elev=%s) alloc(%s perc=%.3f bg=%.3f elev=%.3f) noise(%s perc=%.3f)",
                mActiveSampleCount, fullyAllocated, mFinalBucketGroups.size(), sizeToStr(mTotalCount), sizeToStr(initBackgroundCount),
                sizeToStr(mElevatedCount), sizeToStr(mTotalAllocatedCount), mTotalAllocatedCount/mTotalCount,
                backgroundAlloc/mTotalCount, elevatedAllocatedCount/mTotalCount,
                sizeToStr(noiseAllocated), mElevatedCount > 0 ? noiseAllocated/mElevatedCount : 0));
    }

    public void postRunAnalysis()
    {
        tagMajorGroups();
        logBucketGroups();
        logSampleResults();
        analyseSampleUnallocCounts();
        logWorstAllocatedSamples();
        logSimilarSampleContribs();
        markSimilarGroups();
        logSigReconstructions();
    }

    public void logSampleResults()
    {
        final List<String> categories = mExtSampleData.getFieldNames();

        int fullyAllocCount = 0;
        int partiallyAllocCount = 0;
        int tinyMatchCount = 0;
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
            else if(partiallyAllocated)
                ++partiallyAllocCount;
            else
                ++tinyMatchCount;

            LOGGER.debug(String.format("sample(%d: %s) %s allocated: groups(%d) buckets(%d unalloc=%d) count(total=%s bg=%s elev=%s alloc=%.3f noise=%.2f of %s) cancer(%s) effects(%d: %s)",
                    sampleId, sample.name(), fullyAllocated ? "fully" : (partiallyAllocated ? "partially" : "tiny"),
                    sample.getBucketGroups().size(), samBucketList.size(), sample.getUnallocBuckets().size(),
                    sizeToStr(sample.getTotalCount()), sizeToStr(bgTotal),
                    sizeToStr(sample.getElevatedCount()), sample.getAllocPercent(), sample.getNoisePerc(), sizeToStr(sample.getMaxNoise()),
                    cancerType, effectsCount, effects));
        }

        double percAllocated = countAllocated / mElevatedCount;

        LOGGER.debug(String.format("sample summary: total(%d) alloc(%d, %.3f of %s) partial(%d) tiny(%d)",
                mActiveSampleCount, fullyAllocCount, percAllocated, sizeToStr(mElevatedCount),
                partiallyAllocCount, tinyMatchCount));

        logOverallStats();
    }

    public void logBucketGroups()
    {
        List<BucketGroup> bgList = mFinalBucketGroups;

        // sort by score before logging - for now log in order of creation (by highest allocation)
        // Collections.sort(bgList);

        // log top groups
        int maxToLog = bgList.size();

        LOGGER.debug("logging all bucket groups of total({})", bgList.size());

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

            if(bucketGroup.getGroupType().isEmpty())
                linkData = mAnalyser.isMajorGroup(bucketGroup) ? BG_TYPE_MAJOR : BG_TYPE_MINOR;
            else
                linkData = bucketGroup.getGroupType();

            if(!bucketGroup.getTag().isEmpty())
            {
                linkData += String.format(" tag=%s", bucketGroup.getTag());
            }

            String similarGroup = "";
            if(bucketGroup.getMaxSimilarGroup() != null)
            {
                similarGroup = String.format(" simGroup(%d %.3f)", bucketGroup.getMaxSimilarGroup().getId(), bucketGroup.getMaxSimilarScore());
            }

            LOGGER.debug(String.format("rank %d: bg(%d) %s cancer(%s) samples(%d) variants(avg=%s avgAllocPerc=%.3f total=%s perc=%.3f) buckets(%d: %s) effects(%s)%s",
                    i, bucketGroup.getId(), linkData, bucketGroup.getCancerType(), bucketGroup.getSampleIds().size(),
                    sizeToStr(bucketGroup.getAvgCount()), avgAllocPerc, sizeToStr(bucketGroup.getTotalCount()), groupPerc,
                    bucketGroup.getBucketIds().size(), bucketIdsStr, bucketGroup.getEffects(), similarGroup));
        }
    }

    private void analyseSampleUnallocCounts()
    {
        LOGGER.debug("analysing unallocated sample counts");

        for(final SampleData sample : mSampleData)
        {
            if (sample.isExcluded())
                continue;

            if (sample.getAllocPercent() > SAMPLE_ALLOCATED_PERCENT)
                continue;

            final List<BucketGroup> sampleGroupList = sample.getElevBucketGroups();

            final double[] unallocCounts = sample.getUnallocBucketCounts();
            double sampleCount = sample.getElevatedCount();

            double maxCss = 0;
            double maxPotAlloc = 0;
            BucketGroup maxCssGroup = null;
            BucketGroup maxPotAllocGroup = null;

            for (final BucketGroup bucketGroup : mFinalBucketGroups)
            {
                if (bucketGroup.isBackground() || sampleGroupList.contains(bucketGroup))
                    continue;

                final double[] bucketRatios = bucketGroup.getBucketRatios();
                final double[] ratioRanges = bucketGroup.getRatioRanges();

                double[] allocCounts = new double[mBucketCount];
                double allocTotal = sample.getPotentialUnallocCounts(bucketRatios, bucketGroup.getBucketIds(), ratioRanges, allocCounts);

                if (allocTotal > maxPotAlloc)
                {
                    maxPotAlloc = allocTotal;
                    maxPotAllocGroup = bucketGroup;
                }

                double css = calcCosineSim(unallocCounts, bucketGroup.getBucketRatios());

                if (css > maxCss)
                {
                    maxCss = css;
                    maxCssGroup = bucketGroup;
                }
            }

            if (maxPotAlloc / sampleCount >= MIN_GROUP_ALLOC_PERCENT * 1.05
            && maxPotAlloc >= mAnalyser.getMinSampleAllocCount() && maxPotAllocGroup != null)
            {
                LOGGER.debug(String.format("sample(%d) current alloc(%.3f of %s) missed maxAlloc(%s %.3f) with bg(%d)",
                        sample.Id, sample.getAllocPercent(), sizeToStr(sampleCount), sizeToStr(maxPotAlloc),
                        maxPotAlloc / sampleCount, maxPotAllocGroup.getId()));
            }

            if (sample.getAllocPercent() < (1 - MIN_GROUP_ALLOC_PERCENT) && maxCss >= 0.98 && maxCssGroup != null)
            {
                LOGGER.debug(String.format("sample(%d) current alloc(%.3f of %s) high maxCss(%.3f) with to bg(%d)",
                        sample.Id, sample.getAllocPercent(), sizeToStr(sampleCount), maxCss, maxCssGroup.getId()));
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
                    worstIndex, sample.Id, sample.name(), sample.cancerType(), sizeToStr(unallocTotal),
                    sizeToStr(sample.getElevatedCount()), sample.getUnallocPercent(), unallocTotal/mElevatedCount, sample.getBucketGroups().size()));
        }

        double allUnallocTotal = mTotalCount - mTotalAllocatedCount;

        LOGGER.debug(String.format("worst %d (perc=%.3f) samples have unalloc(%s) ofAllUnalloc(perc=%.3f of %s) ofTotal(%.3f)",
                sampleCount, sampleCount/(double)nonExcludedSampleCount, sizeToStr(totalUnallocated), totalUnallocated/allUnallocTotal,
                sizeToStr(allUnallocTotal), totalUnallocated/mTotalAllocatedCount));

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

        Matrix contribMatrix = new Matrix(sigCount, mSampleCount);
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

        // exclude sigs with too few samples
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

                List<Integer> sam1BgList = sample1.getElevBucketGroups().stream().filter(x -> x.isBackground() == false).map(BucketGroup::getId).collect(Collectors.toList());
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

    public void tagMajorGroups()
    {
        for (final BucketGroup bucketGroup : mFinalBucketGroups)
        {
            if(bucketGroup.isBackground())
                continue;

            if(!bucketGroup.getGroupType().isEmpty())
                continue;

            if(mAnalyser.isMajorGroup(bucketGroup))
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

        // LOGGER.debug("testing sig reconstruction");

        // test whether any sigs can be reconstructed from the others
        CountsSigContribOptimiser sigOptim = new CountsSigContribOptimiser(mBucketCount, false, 0.99);

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

            if (sigOptim.getAllocPerc() < 0.9)
            {
                if(sigOptim.getAllocPerc() >= 0.5)
                {
                    LOGGER.debug(String.format("bg(%d) achieved low reconstruction from %d sigs to %.3f percent",
                            testGroup.getId(), sigOptim.contributingSigCount(), sigOptim.getAllocPerc()));
                }

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

                    if(otherGroup.getMaxSimilarGroup() != testGroup) // already known
                        otherGroup.addGroupLinks(String.format("corr_bg_%d", testGroup.getId()));
                }
            }
        }
    }

    private void markSimilarGroups()
    {
        // add annotations for very similar group
        for (BucketGroup bucketGroup : mFinalBucketGroups)
        {
            if (bucketGroup.getMaxSimilarGroup() == null)
                continue;

            if (bucketGroup.getMaxSimilarScore() < SIG_SIMILAR_CSS * 0.9)
                continue;

            bucketGroup.addGroupLinks(String.format("sim_bg_%d_%.2f",
                    bucketGroup.getMaxSimilarGroup().getId(), bucketGroup.getMaxSimilarScore()));
        }
    }

    public void compareSignatures()
    {
        double internalSigCssThreshold = SIG_SIMILAR_CSS * 0.9;

        if(mProposedSigs != null && !mUsingRefSigs)
        {
            // first the internally generated ones
            List<double[]> cssResults = getTopCssPairs(mProposedSigs, mProposedSigs, internalSigCssThreshold, true, true, true, false);

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
                    if(bg1.isBackground() || bg2.isBackground())
                        continue;

                    LOGGER.debug(String.format("proposed sig(%s bg=%d: ct=%s eff=%s samples=%d) matches sig(%d bg=%d: ct=%s eff=%s samples=%d) with css(%.4f)",
                            sigId1, bg1.getId(), bg1.getCancerType(), bg1.getEffects(), bg1.getSampleIds().size(),
                            sigId2, bg2.getId(), bg2.getCancerType(), bg2.getEffects(), bg2.getSampleIds().size(), css));

                    // skip these annotation if already known
                    if(bg1.getMaxSimilarGroup() != bg2)
                        bg1.addGroupLinks(String.format("css_bg_%d_%.2f", bg2.getId(), css));

                    if(bg2.getMaxSimilarGroup() != bg1)
                        bg2.addGroupLinks(String.format("css_bg_%d_%.2f", bg1.getId(), css));
                }
            }
        }

        if(mReferenceSigs == null || mProposedSigs == null)
            return;

        double externalSigCssThreshold = 0.9;

        List<double[]> cssResults = getTopCssPairs(mReferenceSigs, mProposedSigs, externalSigCssThreshold, false, false);

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

}
