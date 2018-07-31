package com.hartwig.hmftools.data_analyser.calcs;

import static java.lang.Math.ceil;
import static java.lang.Math.floor;
import static java.lang.Math.log10;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;

import static com.hartwig.hmftools.data_analyser.DataAnalyser.OUTPUT_DIR;
import static com.hartwig.hmftools.data_analyser.DataAnalyser.OUTPUT_FILE_ID;
import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.CSSR_I1;
import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.CSSR_I2;
import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.CSSR_VAL;
import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.calcCSS;
import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.getLeastSimilarEntry;
import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.getTopCssPairs;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getNewFile;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.listToArray;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sumVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.writeMatrixData;
import static com.hartwig.hmftools.data_analyser.calcs.NmfConfig.NMF_REF_SIG_FILE;
import static com.hartwig.hmftools.data_analyser.types.BucketGroup.getCombinedBuckets;
import static com.hartwig.hmftools.data_analyser.types.BucketGroup.getMatchingBucketList;
import static com.hartwig.hmftools.data_analyser.types.GenericDataCollection.GD_TYPE_STRING;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.PerformanceCounter;
import com.hartwig.hmftools.data_analyser.loaders.GenericDataLoader;
import com.hartwig.hmftools.data_analyser.types.BucketGroup;
import com.hartwig.hmftools.data_analyser.types.GenericDataCollection;
import com.hartwig.hmftools.data_analyser.types.NmfMatrix;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class BucketAnalyser {

    private static final Logger LOGGER = LogManager.getLogger(NmfManager.class);

    private GenericDataCollection mDataCollection;

    private NmfMatrix mSampleCounts;

    private double[] mSampleTotals;
    private NmfMatrix mBucketMedianRatios;
    private NmfMatrix mBucketProbs;
    private NmfMatrix mBackgroundCounts;
    HashMap<Integer, List<Integer>> mAllSampleBucketGroups;
    HashMap<Integer, List<Integer>> mWorkingSBGroups; // will be reduced as samples are assigned
    private NmfMatrix mProposedSigs;
    private NmfMatrix mReferenceSigs;

    // external sample data for verification and correlation
    GenericDataCollection mExtSampleData;
    HashMap<String,Integer> mExtCategoriesMap;

    private List<BucketGroup> mBucketGroups;

    // config
    private String mOutputDir;
    private String mOutputFileId;
    private double mHighCssThreshold; // CSS level for samples or groups to be consider similar
    private double mLowerCssThreshold; // applied once clear bucket groups have been created
    private double mHighRequiredMatch; // high-level required number of buckets to match between samples or groups
    private double mLowerRequiredMatch; // lower-level match
    private int mMaxProposedSigs;

    private static String BA_CSS_HIGH_THRESHOLD = "ba_css_high";
    private static String BA_CSS_LOW_THRESHOLD = "ba_css_low";
    private static String BA_MAX_PROPOSED_SIGS = "ba_max_proposed_sigs";

    private static int MIN_BG_BUCKET_COUNT = 20; // below which a bucket ratio is considered invalid
    private static double MIN_BUCKET_ALLOCATION = 0.9; // consider a sample fully allocated if X% of its buckets are attributed to groups
    private static double DOMINANT_CATEGORY_PERCENT = 0.75; /// mark a group with a category if X% of samples in it have this attribute (eg cancer type, UV)

    // external data file attributes
    private static String BA_EXT_SAMPLE_DATA_FILE = "ba_ext_data_file";
    private static int COL_SAMPLE_ID = 0;
    private static int COL_CANCER_TYPE = 1;
    private static int CATEGORY_COL_COUNT = 3;
    private static String CATEGORY_CANCER_TYPE = "Cancer";
    private static String EFFECT_TRUE = "TRUE";

    private List<Integer> mSampleWatchList;

    public BucketAnalyser()
    {
        mOutputDir = "";
        mOutputFileId = "";

        mDataCollection = null;
        mSampleCounts = null;
        mBucketMedianRatios = null;
        mSampleTotals = null;
        mBackgroundCounts = null;
        mExtSampleData = null;
        mExtCategoriesMap = null;
        mReferenceSigs = null;
        mProposedSigs = null;

        mBucketGroups = Lists.newArrayList();

        mHighCssThreshold = 0;
        mLowerCssThreshold = 0;
        mHighRequiredMatch = 0.95;
        mLowerRequiredMatch = 0.3;
        mMaxProposedSigs = 0;

        mSampleWatchList = Lists.newArrayList();

        //mSampleWatchList.add(1424);
    }

    public static void addCmdLineArgs(Options options) {

        options.addOption(BA_EXT_SAMPLE_DATA_FILE, true, "Sample external data");
        options.addOption(BA_CSS_HIGH_THRESHOLD, true, "Cosine sim for high-match test");
        options.addOption(BA_CSS_LOW_THRESHOLD, true, "Cosine sim for low-match test");
        options.addOption(BA_MAX_PROPOSED_SIGS, true, "Maximum number of bucket groups to turn into proposed sigs");
    }

    public void initialise(GenericDataCollection collection, final CommandLine cmd)
    {
        mDataCollection = collection;
        mOutputFileId = cmd.getOptionValue(OUTPUT_FILE_ID);
        mOutputDir = cmd.getOptionValue(OUTPUT_DIR);

        mHighCssThreshold = Double.parseDouble(cmd.getOptionValue(BA_CSS_HIGH_THRESHOLD, "0.995"));
        mLowerCssThreshold = Double.parseDouble(cmd.getOptionValue(BA_CSS_LOW_THRESHOLD, "0.98"));
        mHighRequiredMatch = 0.95;
        mLowerRequiredMatch = 0.3;
        mMaxProposedSigs = Integer.parseInt(cmd.getOptionValue(BA_MAX_PROPOSED_SIGS, "0"));

        LOGGER.info("config: cssThreshold(high={} low={} reqMatch(high={} low={}",
                mHighCssThreshold, mLowerCssThreshold, mHighRequiredMatch, mLowerRequiredMatch);

        mSampleCounts = DataUtils.createMatrixFromListData(mDataCollection.getData());
        mSampleCounts.cacheTranspose();

        if(cmd.hasOption(NMF_REF_SIG_FILE))
        {
            GenericDataCollection dataCollection = GenericDataLoader.loadFile(cmd.getOptionValue(NMF_REF_SIG_FILE));
            mReferenceSigs = DataUtils.createMatrixFromListData(dataCollection.getData());
            mReferenceSigs.cacheTranspose();
        }

        LOGGER.info("bucketCount({}) sampleCount({})", mSampleCounts.Rows, mSampleCounts.Cols);

        // report bucket info
        for (int i = 0; i < mSampleCounts.Rows; ++i)
        {
            int bucketTotal = (int) DataUtils.sumVector(mSampleCounts.getRow(i));
            // LOGGER.debug("bucket({}) count({})", i, bucketTotal);
        }

        if(cmd.hasOption(BA_EXT_SAMPLE_DATA_FILE))
        {
            final String fileName = cmd.getOptionValue(BA_EXT_SAMPLE_DATA_FILE);
            mExtSampleData = GenericDataLoader.loadFile(fileName, GD_TYPE_STRING);

            List<Integer> sampleIds = Lists.newArrayList();
            for (int i = 0; i < mSampleCounts.Cols; ++i)
            {
                sampleIds.add(i);
            }

            mExtCategoriesMap = populateSampleCategoryMap(sampleIds);

            LOGGER.debug("registered {} cancer types and sub-categories", mExtCategoriesMap.size());
        }
    }

    public void run()
    {
        PerformanceCounter perfCounter = new PerformanceCounter("BucketMeanRatios");

        perfCounter.start("CalcRatios");
        calcBucketMeanRatios();
        collectElevatedSampleBuckets();
        perfCounter.stop();

        perfCounter.start("BucketGroups");

        // start by looking for near-exact matches in buckets across samples, forming bucket groups
        formBucketGroups(mHighRequiredMatch);

        broadenBucketGroupDefinitions();

        checkBucketsAgainstAllSamples();

        // in this version, a sample can be part of more than one group with overlapping buckets
        // so its counts can be included more than one, so some latter recalc will be required
        // befor sigs are inferred from the bucket count ratios
        mergeBucketGroups(mHighCssThreshold);

        /*
        // now consider samples belonging to more than 1 bucket group
        multiAssignSamples(mLowerRequiredMatch);

        // this time lower the common bucket counts in an attempt to incease the CSS
        formSampleSubGroups();

        mergeBucketGroups(mLowerCssThreshold);

        // try again - is anything more added?
        multiAssignSamples(mLowerRequiredMatch);
        */

        perfCounter.stop();

        logBucketGroups(false);

        perfCounter.start("AnalyseResults");
        analyseGroupsVsExtData();
        logBucketGroups(true);

        logSampleResults();

        createSignatures();
        compareSignatures();
        perfCounter.stop();

        writeSampleData();

        // split each sample's counts into background and deregulated counts
        perfCounter.start("SplitCounts");
        calcBackgroundSampleCounts();
        splitSampleCounts();
        perfCounter.stop();

        perfCounter.logStats();
    }

    private void calcBucketMeanRatios()
    {
        int bucketCount = mSampleCounts.Rows;
        int sampleCount = mSampleCounts.Cols;

        mSampleTotals = new double[sampleCount];

        for (int i = 0; i < mSampleCounts.Cols; ++i)
        {
            mSampleTotals[i] = sumVector(mSampleCounts.getCol(i));
        }

        // work out bucket median values (literally 50th percentile values
        mBucketMedianRatios = new NmfMatrix(bucketCount, sampleCount);
        mBucketProbs = new NmfMatrix(bucketCount, sampleCount);

        double[] medianBucketCounts = new double[bucketCount];

        double bucketMedRange = 0.005;
        int minMedIndex = (int) floor(sampleCount * (0.5 - bucketMedRange));
        int maxMedIndex = (int) ceil(sampleCount * (0.5 + bucketMedRange));
        int gridSize = sampleCount * bucketCount;

        LOGGER.debug("calculating median counts from indices({} -> {})", minMedIndex, maxMedIndex);

        for (int i = 0; i < bucketCount; ++i)
        {
            final double[] bucketCounts = mSampleCounts.getRow(i);

            final List<Integer> bcSorted = DataUtils.getSortedVectorIndices(bucketCounts, true);

            double countTotal = 0;
            for (int j = minMedIndex; j <= maxMedIndex; ++j)
            {
                countTotal += bucketCounts[bcSorted.get(j)];
            }

            double medBucketCount = countTotal / (maxMedIndex - minMedIndex + 1);

            LOGGER.debug(String.format("bucket(%d) median count(%.0f)", i, medBucketCount));

            medianBucketCounts[i] = medBucketCount;
        }

        double totalMedCount = sumVector(medianBucketCounts);

        double[][] brData = mBucketMedianRatios.getData();
        double[][] probData = mBucketProbs.getData();
        final double[][] scData = mSampleCounts.getData();

        int probExpMax = 15;
        int[] probFrequencies = new int[probExpMax+1];
        double minProb = pow(10, -probExpMax);
        int zeroProbIndex = probFrequencies.length-1;

        for (int i = 0; i < bucketCount; ++i)
        {
            // percent of this bucket vs total in median terms
            double bucketMedianRatio = medianBucketCounts[i] / totalMedCount;

            for (int j = 0; j < sampleCount; ++j)
            {
                double expMedianCount = bucketMedianRatio * mSampleTotals[j];
                int sbCount = (int)scData[i][j];
                brData[i][j] = sbCount / expMedianCount;

                double prob = 1;

                if(sbCount > expMedianCount)
                {
                    PoissonDistribution poisDist = new PoissonDistribution(expMedianCount);
                    prob = 1 - poisDist.cumulativeProbability(sbCount - 1);
                    prob = min(prob * gridSize, 1);
                }

                probData[i][j] = prob;

                if(prob > minProb && prob < 0.1)
                {
                    int baseProb = -(int)round(log10(prob));

                    if (baseProb >= 0 && baseProb <= probExpMax)
                        probFrequencies[baseProb] += 1;
                }
                else if(prob < minProb)
                {
                    // allocate to the last slot
                    probFrequencies[zeroProbIndex] += 1;
                }
            }
        }

        for (int i = 0; i < probFrequencies.length-1; ++i)
        {
            LOGGER.debug(String.format("probability(1e-%d) freq(%d) percOfTotal(%.4f)",
                    i, probFrequencies[i], probFrequencies[i]/(double)gridSize));
        }

        LOGGER.debug(String.format("probability(zero) freq(%d) percOfTotal(%.4f)",
                probFrequencies[zeroProbIndex], probFrequencies[zeroProbIndex]/(double)gridSize));
    }

    private static double MAX_ELEVATED_PROB = 1e-12;

    private void collectElevatedSampleBuckets()
    {
        int bucketCount = mSampleCounts.Rows;
        int sampleCount = mSampleCounts.Cols;

        mAllSampleBucketGroups = new HashMap();
        mWorkingSBGroups = new HashMap();

        int totalCount = 0;
        double[][] probData = mBucketProbs.getData();
        // double[][] brData = mBucketMedianRatios.getData();

        for (int i = 0; i < sampleCount; ++i)
        {
            List<Integer> bucketList = Lists.newArrayList();

            for (int j = 0; j < bucketCount; ++j)
            {
                // if (brData[j][i] < MIN_ELEVATED_RATIO) // previous ratio-based test

                if (probData[j][i] > MAX_ELEVATED_PROB)
                {
                    continue;
                }

                bucketList.add(j);
                ++totalCount;
            }

            if (bucketList.isEmpty())
            {
                continue;
            }

            mAllSampleBucketGroups.put(i, bucketList);
            mWorkingSBGroups.put(i, bucketList);
            // LOGGER.debug("sample({}) has {} elevated buckets", i, bucketList.size());
        }

        LOGGER.debug(String.format("samples with elevated buckets: count(%d perc=%.2f), buckets(%d perc=%.3f)",
                mAllSampleBucketGroups.size(), mAllSampleBucketGroups.size()/(double)sampleCount,
                totalCount, totalCount/(double)(bucketCount*sampleCount)));
    }

    private void formBucketGroups(double reqMatchPercent)
    {
        // forms groups out of samples with similarly elevated buckets
        // returns true if groups were created or added to
        int sampleCount = mSampleCounts.Cols;

        int groupsAdjusted = 0;
        int groupsCreated = 0;

        for (int samIndex1 = 0; samIndex1 < sampleCount; ++samIndex1)
        {
            final List<Integer> bl1 = mWorkingSBGroups.get(samIndex1);

            if (bl1 == null)
                continue;

            boolean addedToGroup = false;

            for (BucketGroup bucketGroup : mBucketGroups)
            {
                final List<Integer> commonBuckets = getMatchingBucketList(bucketGroup.getBucketIds(), bl1);

                double groupMatch = commonBuckets.size() / (double)bucketGroup.getBucketIds().size();
                double sampleMatch = commonBuckets.size() / (double)bl1.size();
                double minMatch = min(groupMatch, sampleMatch);

                if (minMatch >= reqMatchPercent)
                {
                    List<Integer> combinedBuckets = getCombinedBuckets(bl1, bucketGroup.getBucketIds());

                    double[] sc1 = extractBucketCountSubset(samIndex1, combinedBuckets);
                    double bcCss = calcSharedCSS(sc1, bucketGroup.getBucketCounts());

                    if(bcCss >= mHighCssThreshold)
                    {
                        // group's buckets remain the same, neither increased nor refined, just add this new sample's counts
                        bucketGroup.addSample(samIndex1, sc1);
                        ++groupsAdjusted;

                        LOGGER.debug(String.format("bg(%d) added sample(%d) with buckets(grp=%d sam=%d match=%d) css(%.4f) totalSamples(%d)",
                                bucketGroup.getId(), samIndex1, bucketGroup.getBucketIds().size(), bl1.size(),
                                commonBuckets.size(), bcCss, bucketGroup.getSampleIds().size()));

                        addedToGroup = true;
                        break;
                    }
                }
            }

            if (addedToGroup)
            {
                mWorkingSBGroups.remove(samIndex1);
                continue;
            }

            for (int samIndex2 = samIndex1 + 1; samIndex2 < sampleCount; ++samIndex2)
            {
                final List<Integer> bl2 = mWorkingSBGroups.get(samIndex2);

                if (bl2 == null)
                    continue;

                final List<Integer> commonBuckets = getMatchingBucketList(bl1, bl2);

                double sample1Match = commonBuckets.size() / (double)bl1.size();
                double sample2Match = commonBuckets.size() / (double)bl2.size();
                double minMatch = min(sample1Match, sample2Match);

                if (minMatch >= reqMatchPercent)
                {
                    List<Integer> combinedBuckets = getCombinedBuckets(bl1, bl2);
                    double[] sc1 = extractBucketCountSubset(samIndex1, combinedBuckets);
                    double[] sc2 = extractBucketCountSubset(samIndex2, combinedBuckets);
                    double bcCss = calcSharedCSS(sc1, sc2);

                    if(bcCss >= mHighCssThreshold)
                    {
                        BucketGroup bucketGroup = new BucketGroup(mBucketGroups.size());
                        bucketGroup.addBuckets(commonBuckets);
                        bucketGroup.addSample(samIndex1, sc1);
                        bucketGroup.addSample(samIndex2, sc2);

                        LOGGER.debug(String.format("added bg(%d) samples(%d and %d) with buckets(s1=%d and s2=%d matched=%d) css(%.4f)",
                                bucketGroup.getId(), samIndex1, samIndex2, bl1.size(), bl2.size(), commonBuckets.size(), bcCss));

                        mBucketGroups.add(bucketGroup);
                        ++groupsCreated;

                        mWorkingSBGroups.remove(samIndex2);
                        mWorkingSBGroups.remove(samIndex1);
                        break;
                    }
                }
            }
        }

        if(mBucketGroups.isEmpty())
        {
            LOGGER.debug("no bucket groups created");
            return;
        }

        LOGGER.debug("bucket groups created({}) additions({}) total({}), unallocated elevated samples({})",
                groupsCreated, groupsAdjusted, mBucketGroups.size(), mWorkingSBGroups.size());
    }

    private void broadenBucketGroupDefinitions()
    {
        int bucketCount = mSampleCounts.Rows;

        for (BucketGroup bucketGroup : mBucketGroups)
        {
            // starting with the current set of buckets, attempt to add in a new bucket as long as the CSS remains above the threshold
            // but only apply if the samples in the group show evidence of being elevated in the new buckets being considered
            final List<Integer> startBuckets = bucketGroup.getBucketIds();
            final List<Integer> sampleIds = bucketGroup.getSampleIds();
            List<Integer> workingBuckets = bucketGroup.getBucketIds();

            // populate the matrix
            NmfMatrix sampleCounts = new NmfMatrix(bucketCount, sampleIds.size());
            double[][] scData = sampleCounts.getData();
            final double[][] refData = mSampleCounts.getData();
            final double[][] bmrData = mBucketMedianRatios.getData();

            for(int i = 0; i < bucketCount; ++i)
            {
                if(!startBuckets.contains(i))
                    continue;

                int sampleIndex = 0;
                for (Integer sampleId : sampleIds)
                {
                    scData[i][sampleIndex] = refData[i][sampleId];
                    ++sampleIndex;
                }
            }

            List<double[]> cssResults = getTopCssPairs(sampleCounts, sampleCounts, mHighCssThreshold, false, true);
            int bgSampleCount = sampleIds.size();
            int expectedPairs = bgSampleCount * (bgSampleCount-1) / 2;
            int bucketsAdded = 0;

            int lastTestedBucket = -1;
            for (int testBucket = 0; testBucket < bucketCount; ++testBucket)
            {
                if (workingBuckets.contains(testBucket))
                    continue;

                // check if this new bucket tested shows evidence of being elevated in the samples
                boolean areElevated = true;
                for (Integer sampleId : sampleIds)
                {
                    if(bmrData[testBucket][sampleId] < 1)
                    {
                        areElevated = false;
                        break;
                    }
                }

                if(!areElevated)
                    continue;

                // otherwise add this bucket and test out CSS across all samples
                int sampleIndex = 0;
                for (Integer sampleId : sampleIds)
                {
                    // set count for bucket to be tested
                    scData[testBucket][sampleIndex] = refData[testBucket][sampleId];

                    // reset previous
                    if(lastTestedBucket != -1)
                        scData[lastTestedBucket][sampleIndex] = 0;

                    ++sampleIndex;
                }

                lastTestedBucket = testBucket;
                cssResults = getTopCssPairs(sampleCounts, sampleCounts, mHighCssThreshold, false, true);

                if(expectedPairs == cssResults.size())
                {
//                    LOGGER.debug(String.format("bg(%d) added new bucket(%d) worstCss(%.4f) in %d samples",
//                            bucketGroup.getId(), testBucket, cssResults.get(expectedPairs-1)[CSSR_VAL], bgSampleCount));

                    workingBuckets.add(testBucket);
                    ++bucketsAdded;
                    lastTestedBucket = -1; // to retain this for future tests

                    // add this bucket
                    bucketGroup.addBucket(testBucket, sampleCounts.getRow(testBucket));
                }
            }

            if(bucketsAdded > 0)
            {
                LOGGER.debug("bg({}) add {} buckets with {} samples", bucketGroup.getId(), bucketsAdded, bgSampleCount);
            }
        }
    }

    private void checkBucketsAgainstAllSamples()
    {
        LOGGER.debug("checking all samples against bucket groups");

        // take the top X samples and check them against all samples, regardless of whether deregulated
        int samplesAdded = 0;
        int sampleCount = mSampleCounts.Cols;
        int minElevatedBucketMatch = 1;

        final double[][] bmrData = mBucketMedianRatios.getData();

        for (BucketGroup bucketGroup : mBucketGroups)
        {
            if(bucketGroup.calcScore() < 10)
                continue;

            if(bucketGroup.getBucketIds().size() == 1)
                continue; // can't calc CSS using a single bucket

            final List<Integer> groupBuckets = bucketGroup.getBucketIds();

            for (int sampleId = 0; sampleId < sampleCount; ++sampleId)
            {
                if(bucketGroup.hasSample(sampleId))
                    continue; // ignore those already assigned

                double[] samCounts = extractBucketCountSubset(sampleId, groupBuckets);
                double bcCss = calcSharedCSS(samCounts, bucketGroup.getBucketCounts());

                int samElevatedCount = 0;
                int commonBucketCount = 0;
                final List<Integer> samBuckets = mAllSampleBucketGroups.get(sampleId);

                if(samBuckets != null)
                {
                    samElevatedCount = samBuckets != null ? samBuckets.size() : 0;
                    final List<Integer> commonBuckets = getMatchingBucketList(samBuckets, groupBuckets);
                    commonBucketCount = commonBuckets.size();
                }

                if (bcCss >= mHighCssThreshold && commonBucketCount >= minElevatedBucketMatch)
                {
                    bucketGroup.addSample(sampleId, samCounts);
                    ++samplesAdded;

                    // report avg BMR for informational purposes
                    double bmrTotal = 0;
                    for(Integer bucket : groupBuckets)
                    {
                        bmrTotal += bmrData[bucket][sampleId];
                    }

                    double bmrAvg = bmrTotal / groupBuckets.size();

                    LOGGER.debug(String.format("bg(%d) adding sample(%d) with buckets(grp=%d sam=%d match=%d bmrAvg=%.1f) css(%.4f) bg samples(%d)",
                            bucketGroup.getId(), sampleId, groupBuckets.size(), samElevatedCount, commonBucketCount,
                            bmrAvg, bcCss, bucketGroup.getSampleIds().size()));
                }
            }
        }

        LOGGER.debug("added {} samples to existing bucket groups", samplesAdded);
    }


    private int MIN_BUCKET_COUNT_OVERLAP = 3;

    private void formSampleSubGroups()
    {
        // forms groups out of samples with similarly elevated buckets
        // returns true if groups were created or added to
        int sampleCount = mSampleCounts.Cols;

        int groupsAdjusted = 0;
        int groupsCreated = 0;

        for (int samIndex1 = 0; samIndex1 < sampleCount; ++samIndex1)
        {
            final List<Integer> bl1 = mWorkingSBGroups.get(samIndex1);

            if (bl1 == null)
                continue;

            if(mSampleWatchList.contains(samIndex1))
            {
                LOGGER.debug("specific sample");
            }

            boolean sam1Allocated = false;

            for (int samIndex2 = samIndex1 + 1; samIndex2 < sampleCount; ++samIndex2)
            {
                final List<Integer> bl2 = mWorkingSBGroups.get(samIndex2);

                if (bl2 == null)
                    continue;

                final List<Integer> commonBuckets = getMatchingBucketList(bl1, bl2);
                int commonBucketCount = commonBuckets.size();

                if(areAnySampleBucketsAllocated(samIndex1, commonBuckets))
                {
                    // LOGGER.debug("bg({}) vs sample({}) when already allocated to another BG", bucketGroup.getId(), samIndex1);
                    continue;
                }

                if(commonBucketCount <= MIN_BUCKET_COUNT_OVERLAP)
                    continue;

                boolean groupCreated = false;

                while(commonBuckets.size() > MIN_BUCKET_COUNT_OVERLAP)
                {
                    // test out with fewer and fewer buckets to see if a high CSS can be found
                    double[] sc1 = extractBucketCountSubset(samIndex1, commonBuckets);
                    double[] sc2 = extractBucketCountSubset(samIndex2, commonBuckets);
                    double bcCss = calcSharedCSS(sc1, sc2);

                    if(bcCss >= mHighCssThreshold)
                    {
                        BucketGroup bucketGroup = new BucketGroup(mBucketGroups.size());
                        bucketGroup.addBuckets(commonBuckets);
                        bucketGroup.addSample(samIndex1, sc1);
                        bucketGroup.addSample(samIndex2, sc2);

                        LOGGER.debug(String.format("added bg(%d) samples(%d and %d) with buckets(s1=%d and s2=%d matched=%d orig=%d) css(%.4f)",
                                bucketGroup.getId(), samIndex1, samIndex2, bl1.size(), bl2.size(), commonBuckets.size(), commonBucketCount, bcCss));

                        mBucketGroups.add(bucketGroup);
                        groupCreated = true;
                        break;
                    }
                    else if(commonBuckets.size() == MIN_BUCKET_COUNT_OVERLAP + 1)
                    {
                        break;
                    }

//                    LOGGER.debug(String.format("samples(%d and %d) testing with common buckets(%s) vs orig(%d) lastCss(%.4f)",
//                            samIndex1, samIndex2, commonBuckets.size(), commonBucketCount, bcCss));

                    // remove the worst of the common buckets
                    int nextIndex = getLeastSimilarEntry(sc1, sc2);

                    if(nextIndex < 0)
                        break;

                    for(Integer index : commonBuckets)
                    {
                        if(index == nextIndex)
                        {
                            commonBuckets.remove(index);
                            break;
                        }
                    }

                } // end test of CSS with less common buckets

                if(groupCreated)
                {
                    ++groupsCreated;

                    if(areMostSampleBucketsAllocated(samIndex2, bl2, MIN_BUCKET_ALLOCATION))
                    {
                        LOGGER.debug("sample({}) now sufficiently allocated to bucket(s)", samIndex2);
                        mWorkingSBGroups.remove(samIndex2);
                    }

                    if(areMostSampleBucketsAllocated(samIndex1, bl1, MIN_BUCKET_ALLOCATION))
                    {
                        sam1Allocated = true;
                        break; // nothing more to test against sample 1
                    }
                }

            } // end for each sample 2

            if(sam1Allocated)
            {
                LOGGER.debug("sample({}) now sufficiently allocated to bucket(s)", samIndex1);
                mWorkingSBGroups.remove(samIndex1);
            }

        } // end for each sample 1

        LOGGER.debug("bucket groups created({}) additions({}) total({}), unallocated elevated samples({})",
                groupsCreated, groupsAdjusted, mBucketGroups.size(), mWorkingSBGroups.size());
    }

    private void multiAssignSamples(double reqMatchPercent)
    {
        LOGGER.debug("checking multi-assignments for elevated {} samples and {} bucket groups, requiredMatchPerc({})",
                mWorkingSBGroups.size(), mBucketGroups.size(), reqMatchPercent);

        // check unallocated samples against each other and existing groups, looking for partial overlap
        int sampleCount = mSampleCounts.Cols;

        int groupsAdjusted = 0;
        int groupsCreated = 0;

        for (int samIndex1 = 0; samIndex1 < sampleCount; ++samIndex1)
        {
            if(mSampleWatchList.contains(samIndex1))
            {
                LOGGER.debug("specific sample");
            }

            final List<Integer> bl1 = mWorkingSBGroups.get(samIndex1);

            if (bl1 == null)
                continue;

            boolean addedToGroup = false;

            // work out how well each bucket group covers this sample, and at end assign starting with the best
            double[] assignmentScores = new double[mBucketGroups.size()];
            boolean hasPossibleAssignments = false;

            for (int i = 0; i < mBucketGroups.size(); ++i)
            {
                BucketGroup bucketGroup = mBucketGroups.get(i);

                if(bucketGroup.hasSample(samIndex1))
                    continue;

                final List<Integer> groupBuckets = bucketGroup.getBucketIds();
                final List<Integer> commonBuckets = getMatchingBucketList(groupBuckets, bl1);

                if(areAnySampleBucketsAllocated(samIndex1, commonBuckets))
                {
                    continue;
                }

                /*
                // look for a sample with an excess of elevated buckets
                if (bl1.size() <= commonBuckets.size())
                    continue;

                    // previously checked sample match vs required below, but now allowing samples to be linked to group
                    // even if the group covers more buckets
                */

                // test similarity for the common buckets
                double[] sc1 = extractBucketCountSubset(samIndex1, commonBuckets);
                double bcCss = calcSharedCSS(sc1, bucketGroup.getBucketCounts());

                double sampleMatch = commonBuckets.size() / (double)bl1.size();
                double groupMatch = commonBuckets.size() / (double)groupBuckets.size();
                double minMatch = min(sampleMatch, groupMatch);

                if (minMatch >= reqMatchPercent && bcCss >= mHighCssThreshold)
                {
                    assignmentScores[i] = sampleMatch * bcCss;
                    hasPossibleAssignments = true;
                }
            }

            if(hasPossibleAssignments)
            {
                List<Integer> sortedAssignments = getSortedVectorIndices(assignmentScores, false);

                for(int i = 0; i < sortedAssignments.size(); ++i)
                {
                    int bgIndex = sortedAssignments.get(i);
                    double asgnScore = assignmentScores[bgIndex];

                    if(asgnScore == 0)
                        break;

                    BucketGroup bucketGroup = mBucketGroups.get(bgIndex);

                    final List<Integer> commonBuckets = getMatchingBucketList(bucketGroup.getBucketIds(), bl1);

                    if(addedToGroup && areAnySampleBucketsAllocated(samIndex1, commonBuckets)) // check if has just been assigned
                        continue;

                    double[] sc1 = extractBucketCountSubset(samIndex1, commonBuckets);
                    double bcCss = calcSharedCSS(sc1, bucketGroup.getBucketCounts());

                    // group's buckets remain the same, neither increased nor refined
                    bucketGroup.addSample(samIndex1, sc1);
                    ++groupsAdjusted;
                    addedToGroup = true;

                    LOGGER.debug(String.format("bg(%d) added sample(%d) with buckets(grp=%d sam=%d match=%d) totalSamples(%d) css(%.4f) asgnScore(%.2f)",
                            bucketGroup.getId(), samIndex1, bucketGroup.getBucketIds().size(), bl1.size(),
                            commonBuckets.size(), bucketGroup.getSampleIds().size(), bcCss, asgnScore));
                }
            }

            if (addedToGroup)
            {
                if(areMostSampleBucketsAllocated(samIndex1, bl1, MIN_BUCKET_ALLOCATION))
                {
                    LOGGER.debug("sample({}) now sufficiently allocated to bucket(s)", samIndex1);
                    mWorkingSBGroups.remove(samIndex1);
                    continue;
                }
            }

            // otherwise keep checking amongst unallocated samples in a similar manner
            boolean sam1Allocated = false;

            for (int samIndex2 = samIndex1 + 1; samIndex2 < sampleCount; ++samIndex2)
            {
                final List<Integer> bl2 = mWorkingSBGroups.get(samIndex2);

                if (bl2 == null)
                    continue;

                final List<Integer> commonBuckets = getMatchingBucketList(bl1, bl2);

                if(areAnySampleBucketsAllocated(samIndex1, commonBuckets))
                {
                    continue;
                }

                double[] sc1 = extractBucketCountSubset(samIndex1, commonBuckets);
                double[] sc2 = extractBucketCountSubset(samIndex2, commonBuckets);
                double bcCss = calcSharedCSS(sc1, sc2);

                double sample1Match = commonBuckets.size() / (double)bl1.size();
                double sample2Match = commonBuckets.size() / (double)bl2.size();
                double minMatch = min(sample1Match, sample2Match);

                if(minMatch > 1 && minMatch >= reqMatchPercent && bcCss >= mHighCssThreshold)
                {
                    BucketGroup bucketGroup = new BucketGroup(mBucketGroups.size());
                    bucketGroup.addBuckets(commonBuckets);
                    bucketGroup.addSample(samIndex1, sc1);
                    bucketGroup.addSample(samIndex2, sc2);

                    LOGGER.debug(String.format("added bg(%d) samples(%d and %d) with buckets(s1=%d and s2=%d matched=%d) css(%.4f)",
                            bucketGroup.getId(), samIndex1, samIndex2, bl1.size(), bl2.size(), commonBuckets.size(), bcCss));

                    mBucketGroups.add(bucketGroup);
                    ++groupsCreated;

                    if(areMostSampleBucketsAllocated(samIndex2, bl2, MIN_BUCKET_ALLOCATION))
                    {
                        LOGGER.debug("sample({}) now sufficiently allocated to bucket(s)", samIndex2);
                        mWorkingSBGroups.remove(samIndex2);
                    }

                    if(areMostSampleBucketsAllocated(samIndex1, bl1, MIN_BUCKET_ALLOCATION))
                    {
                        sam1Allocated = true;
                        break;
                    }
                }
            }

            if(sam1Allocated)
            {
                LOGGER.debug("sample({}) now sufficiently allocated to bucket(s)", samIndex1);
                mWorkingSBGroups.remove(samIndex1);
            }
        }

        if(groupsAdjusted == 0 && groupsCreated == 0)
        {
            LOGGER.debug("no new/adjusted bucket groups");
            return;
        }

        LOGGER.debug("bucket groups created({}) adjusted({}) total({}), unallocated elevated samples({})",
                groupsCreated, groupsAdjusted, mBucketGroups.size(), mWorkingSBGroups.size());
    }

    private void mergeBucketGroups(double cssThreshold)
    {
        double reqMatchPercent = 0.85;

        int bgIndex1 = 0;
        while (bgIndex1 < mBucketGroups.size())
        {
            BucketGroup bg1 = mBucketGroups.get(bgIndex1);
            boolean removeBg1 = false;

            int bgIndex2 = bgIndex1+1;
            while (bgIndex2 < mBucketGroups.size())
            {
                BucketGroup bg2 = mBucketGroups.get(bgIndex2);
                boolean removeBg2 = false;

                final List<Integer> commonBuckets = getMatchingBucketList(bg1.getBucketIds(), bg2.getBucketIds());

                int commmonBucketCount = commonBuckets.size();
                double bg1Match = commmonBucketCount / (double) bg1.getBucketIds().size();
                double bg2Match = commmonBucketCount / (double) bg2.getBucketIds().size();
                double minMatch = min(bg1Match, bg2Match);
                double bcCss = calcSharedCSS(bg1.getBucketCounts(), bg2.getBucketCounts());

                if (minMatch >= reqMatchPercent && bcCss >= cssThreshold)
                {
                    int bg1SC = bg1.getSampleIds().size();
                    int bg2SC = bg2.getSampleIds().size();

                    if(minMatch < 1)
                    {
                        // should new buckets be merged in? or just go with the common set or most prevalent set
                        if (bg1SC > 1.5 * bg2SC)
                        {
                            bg1.merge(bg2.getSampleIds(), bg2.getBucketCounts());
                            removeBg2 = true;
                        }
                        else if (bg2SC > 1.5 * bg1SC)
                        {
                            bg2.merge(bg1.getSampleIds(), bg1.getBucketCounts());
                            removeBg1 = true;
                        }
                        else
                        {
                            // go with the common set only
                            bg1.merge(bg2.getSampleIds(), bg2.getBucketCounts());
                            bg1.reduceToBucketSet(commonBuckets);
                            removeBg2 = true;
                        }
                    }
                    else
                    {
                        bg1.merge(bg2.getSampleIds(), bg2.getBucketCounts());
                        removeBg2 = true;
                    }

                    LOGGER.debug(String.format("bg(%d) %s bg(%d) buckets(bg1=%d bg2=%d match=%d) css(%.4f) samples(bg1=%d bg2=%d total=%d)",
                            bg1.getId(), removeBg2 ? "merges in" : "merged into", bg2.getId(), bg1.getBucketIds().size(), bg2.getBucketIds().size(),
                            commmonBucketCount, bcCss, bg1SC, bg2SC, bg1SC + bg2SC));

                    if(removeBg2)
                        mBucketGroups.remove(bgIndex2);
                    else
                        break;

                    continue;
                }
                else if(minMatch >= reqMatchPercent * 0.8 && bcCss >= 0.98)
                {
                    LOGGER.debug(String.format("bg(%d) vs bg(%d) close-to-merge buckets(bg1=%d bg2=%d match=%d) css(%.4f)",
                            bg1.getId(), bg2.getId(), bg1.getBucketIds().size(), bg2.getBucketIds().size(), commmonBucketCount, bcCss));
                }
                else if(minMatch >= 0.3 && bcCss < 0.5)
                {
                    // check for overlapping different processes?
                    LOGGER.debug(String.format("bg(%d) vs bg(%d) very diff despite buckets(bg1=%d bg2=%d match=%d) css(%.4f)",
                            bg1.getId(), bg2.getId(), bg1.getBucketIds().size(), bg2.getBucketIds().size(), commmonBucketCount, bcCss));
                }

                if(removeBg2)
                    mBucketGroups.remove(bgIndex2);
                else
                    ++bgIndex2;

                if(removeBg1)
                    break;
            }

            if(removeBg1)
                mBucketGroups.remove(bgIndex1);
            else
                ++bgIndex1;
        }
    }

    private void logBucketGroups(boolean clearTypesOnly)
    {
        // sort by score before logging
        Collections.sort(mBucketGroups);

        // log top groups
        int maxToLog = clearTypesOnly ? mBucketGroups.size() : min(mBucketGroups.size(), 40);
        int minScore = clearTypesOnly ? 0 : 10; // 4 samples x 5 buckets, or 2 x 10

        LOGGER.debug("logging top {} buckets", clearTypesOnly ? "clear-type" : maxToLog);

        for (int i = 0; i < maxToLog; ++i)
        {
            BucketGroup bucketGroup = mBucketGroups.get(i);

            if (bucketGroup.calcScore() < minScore)
                break;

            if(clearTypesOnly)
            {
                if(bucketGroup.getCancerType().isEmpty())
                    continue;

                LOGGER.debug(String.format("rank %d: bg(%d) cancer(%s) score(%.0f) samples(%d) buckets(%d: %s) effects(%s)",
                        i, bucketGroup.getId(), bucketGroup.getCancerType(), bucketGroup.calcScore(), bucketGroup.getSampleIds().size(),
                        bucketGroup.getBucketIds().size(), bucketGroup.getBucketIds().toString(), bucketGroup.getEffects()));
            }
            else
            {
                LOGGER.debug(String.format("rank %d: bg(%d) score(%.0f) samples(%d) buckets(%d)",
                        i, bucketGroup.getId(), bucketGroup.calcScore(), bucketGroup.getSampleIds().size(),
                        bucketGroup.getBucketIds().size()));
            }
        }
    }

    private void logSampleResults()
    {
        final List<String> categories = mExtSampleData.getFieldNames();

        int fullyAllocated = 0;
        int partiallyAllocated = 0;
        int noMatchCount = 0;
        int someMatchCount = 0;

        for(int sampleId = 0; sampleId < mSampleCounts.Cols; ++sampleId)
        {
            final List<Integer> samBucketList = mAllSampleBucketGroups.get(sampleId);

            if(samBucketList == null)
                continue;

            final String sampleName = getSampleName(sampleId);

            final List<String> sampleData = getSampleExtData(sampleName);
            if(sampleData == null)
                continue;

            if(mSampleWatchList.contains(sampleId))
            {
                LOGGER.debug("specific sample");
            }

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

            int bucketsAllocated = getSampleBucketsAllocated(sampleId, samBucketList);
            int groupCount = getSampleGroupMembershipCount(sampleId);
            double allocPerc = bucketsAllocated/(double)samBucketList.size();

            if(allocPerc >= MIN_BUCKET_ALLOCATION)
            {
                ++fullyAllocated;
                LOGGER.debug(String.format("sample(%d: %s) fully allocated to %d group(s): buckets(%d alloc=%d perc=%.2f) cancer(%s) effects(%d: %s)",
                        sampleId, sampleName, groupCount, samBucketList.size(), bucketsAllocated, allocPerc, cancerType, effectsCount, effects));
                continue;
            }

            if(groupCount > 0)
            {
                ++partiallyAllocated;
                LOGGER.debug(String.format("sample(%d: %s) partially allocated to %d group(s): buckets(%d alloc=%d perc=%.2f) cancer(%s) effects(%d: %s)",
                        sampleId, sampleName, groupCount, samBucketList.size(), bucketsAllocated, allocPerc, cancerType, effectsCount, effects));
                continue;
            }

            // now find the closest matching bucket(s) and try to work out why they didn't match
            boolean someMatch = false;
            for(final BucketGroup bucketGroup : mBucketGroups)
            {
                if(!bucketGroup.getCancerType().equals(cancerType))
                    continue;

                final List<Integer> commonBuckets = getMatchingBucketList(bucketGroup.getBucketIds(), samBucketList);

                if(commonBuckets.size() < 0.5 * samBucketList.size())
                    continue;

                double[] sc1 = extractBucketCountSubset(sampleId, commonBuckets);
                double bcCss = calcSharedCSS(sc1, bucketGroup.getBucketCounts());

                if(bcCss >= 0.9)
                {
                    someMatch = true;

                    LOGGER.debug(String.format("sample(%d: %s ct=%s eff=%d: %s) unallocated close match with bg(%d eff=%s), buckets(grp=%d sam=%d match=%d) css(%.4f)",
                            sampleId, sampleName, cancerType, effectsCount, effects, bucketGroup.getId(), bucketGroup.getEffects(),
                            bucketGroup.getBucketIds().size(), samBucketList.size(), commonBuckets.size(), bcCss));
                }
            }

            if(someMatch)
            {
                ++someMatchCount;
            }
            else
            {
                LOGGER.debug("sample({}: {}) unallocated with no match: buckets({}) cancer({}) effects({}: {})",
                        sampleId, sampleName, samBucketList.size(), cancerType, effectsCount, effects);
                ++noMatchCount;
            }
        }

        int unallocated = someMatchCount + noMatchCount;

        LOGGER.debug("sample summary: total({}) elevated({}) allocated({}) partially({}) unallocated({} possibleMatch={} noMatch={}) groupCount({})",
                mSampleCounts.Cols, mAllSampleBucketGroups.size(), fullyAllocated, partiallyAllocated,
                unallocated, someMatchCount, noMatchCount, mBucketGroups.size());
    }

    private void analyseGroupsVsExtData()
    {
        if(mExtSampleData == null)
            return;

        LOGGER.debug("analysing {} bucket groups", mBucketGroups.size());

        // check counts for each external data category against the samples in each bucket group

        // want to know if:
        // a) the samples have 1 or more identical categoroes and
        // b) the proportion of these categories out of the cohort (ie missing other samples, linked to other bucket groups)

        for (BucketGroup bucketGroup : mBucketGroups)
        {
            final List<Integer> sampleIds = bucketGroup.getSampleIds();
            int sampleCount = sampleIds.size();

            HashMap<String,Integer> categoryCounts = populateSampleCategoryMap(sampleIds);

            for(Map.Entry<String,Integer> entrySet : categoryCounts.entrySet())
            {
                final String catName = entrySet.getKey();
                int catSampleCount = entrySet.getValue();
                double samplesPerc = catSampleCount/(double)sampleCount;

                if(samplesPerc < DOMINANT_CATEGORY_PERCENT * 0.5)
                    continue;

                int catAllCount = mExtCategoriesMap.get(catName);
                double catPerc = catSampleCount / (double) catAllCount;

                LOGGER.debug(String.format("bg(%d) category(%s) count(%d of %d) perc(group=%.3f category=%.3f) buckets(%d)",
                        bucketGroup.getId(), catName, catSampleCount, sampleCount, samplesPerc, catPerc, bucketGroup.getBucketIds().size()));

                if(samplesPerc >= DOMINANT_CATEGORY_PERCENT)
                {
                    if(extractCategoryName(catName).equals(CATEGORY_CANCER_TYPE))
                    {
                        bucketGroup.setCancerType(extractCategoryValue(catName));
                    }
                    else
                    {
                        String effects = bucketGroup.getEffects();

                        if(effects.isEmpty())
                            effects = catName;
                        else
                            effects += ";" + catName;

                        bucketGroup.setEffects(effects);
                    }
                }
            }
        }
    }

    private void compareSignatures()
    {
        if(mProposedSigs != null)
        {
            // first the internally generated ones
            List<double[]> cssResults = getTopCssPairs(mProposedSigs, mProposedSigs, 0.99, true, true);

            if (cssResults.isEmpty())
            {
                LOGGER.debug("no similar proposed sigs from bucket groups");
            } else
            {
                for (final double[] result : cssResults)
                {
                    int bgId1 = (int)result[CSSR_I1];
                    int bgId2 = (int)result[CSSR_I2];
                    final BucketGroup bg1 = mBucketGroups.get(bgId1);
                    final BucketGroup bg2 = mBucketGroups.get(bgId2);

                    LOGGER.debug(String.format("proposed sigs bg(%d: ct=%s eff=%s samples=%d) matches bg(%d: ct=%s eff=%s samples=%d) with css(%.4f)",
                            bg1.getId(), bg1.getCancerType(), bg1.getEffects(), bg1.getSampleIds().size(),
                            bg2.getId(), bg2.getCancerType(), bg2.getEffects(), bg2.getSampleIds().size(), result[CSSR_VAL]));
                }
            }
        }

        if(mReferenceSigs == null)
            return;

        List<double[]> cssResults = getTopCssPairs(mReferenceSigs, mProposedSigs, 0.95, false, false);

        if (cssResults.isEmpty())
        {
            LOGGER.debug("no similar sigs between external ref and bucket groups");
        }
        else
        {
            for (final double[] result : cssResults)
            {
                int externalSigId = (int)result[CSSR_I1] + 1; // bumped up to correspond to convention of starting with 1
                int bgId = (int)result[CSSR_I2];
                final BucketGroup bucketGroup = mBucketGroups.get(bgId);

                LOGGER.debug(String.format("external ref sig(%d) matches bg(%d: ct=%s eff=%s samples=%d) with css(%.4f)",
                        externalSigId, bucketGroup.getId(), bucketGroup.getCancerType(), bucketGroup.getEffects(),
                        bucketGroup.getSampleIds().size(), result[CSSR_VAL]));
            }
        }
    }

    private void calcBackgroundSampleCounts()
    {
        // similar process to elevated counts, but only included sample bucket counts if not elevated
        int bucketCount = mSampleCounts.Rows;
        int sampleCount = mSampleCounts.Cols;

        mBackgroundCounts = new NmfMatrix(bucketCount, sampleCount);
        final double[][] bgData = mBackgroundCounts.getData();

        double[] bgSampleTotals = new double[sampleCount];
        final double[][] scData = mSampleCounts.getData();
        final double[][] probData = mBucketProbs.getData();

        // work out bucket median values (literally 50th percentile values
        // mBucketMedianRatios = new NmfMatrix(bucketCount, sampleCount);
        double[] medianBucketCounts = new double[bucketCount];

        double bucketMedRange = 0.005;

        // first extract counts per bucket but only if not elevated
        for (int i = 0; i < bucketCount; ++i)
        {
            List<Double> counts = Lists.newArrayList();
            for(int j = 0; j < sampleCount; ++j)
            {
                if(probData[i][j] > MAX_ELEVATED_PROB)
                {
                    counts.add(scData[i][j]);
                }
            }

            final double[] bucketCounts = listToArray(counts);
            final List<Integer> bcSorted = DataUtils.getSortedVectorIndices(bucketCounts, true);

            int minMedIndex = (int) floor(counts.size() * (0.5 - bucketMedRange));
            int maxMedIndex = (int) ceil(counts.size() * (0.5 + bucketMedRange));

            double countTotal = 0;
            for (int j = minMedIndex; j <= maxMedIndex; ++j)
            {
                countTotal += bucketCounts[bcSorted.get(j)];
            }

            double medBucketCount = countTotal / (maxMedIndex - minMedIndex + 1);

            // LOGGER.debug(String.format("bucket(%d) median count(%.0f)", i, medBucketCount));

            medianBucketCounts[i] = medBucketCount;
        }

        double totalMedCount = sumVector(medianBucketCounts);

        // repeat to get background sample totals and then expected counts per sample per bucket
        for(int j = 0; j < sampleCount; ++j)
        {
            List<Double> counts = Lists.newArrayList();
            for (int i = 0; i < bucketCount; ++i)
            {
                if (probData[i][j] > MAX_ELEVATED_PROB)
                {
                    counts.add(scData[i][j]);
                }
            }

            final double[] sampleCounts = listToArray(counts);
            bgSampleTotals[j] = sumVector(sampleCounts);

            for (int i = 0; i < bucketCount; ++i)
            {
                // percent of this bucket vs total in median terms
                double bucketMedianRatio = medianBucketCounts[i] / totalMedCount;
                bgData[i][j] = round(bucketMedianRatio * bgSampleTotals[j]);
            }

            // LOGGER.debug(String.format("sample(%d) background count(%.0f) vs actual(%.0f)", j, bgSampleTotals[j], mSampleTotals[j]));
        }
    }

    private void splitSampleCounts()
    {
        LOGGER.debug("splitting sample counts into background and elevated");

        int sampleCount = mSampleCounts.Cols;
        int bucketCount = mSampleCounts.Rows;

        NmfMatrix elevatedData = new NmfMatrix(bucketCount, sampleCount);
        double[][] evData = elevatedData.getData();
        double[][] bgData = mBackgroundCounts.getData();
        double[][] scData = mSampleCounts.getData();

        for(int sampleId = 0; sampleId < sampleCount; ++sampleId)
        {
            final List<Integer> bucketGroups = getSampleBucketGroups(sampleId);
            final List<Integer> elevatedBuckets = mAllSampleBucketGroups.get(sampleId);

            int expAboveActual = 0;
            int elevatedCount = 0;
            int lowProbCount = 0;
            int lowProbExpAboveActual = 0;
            double expAboveActualTotal = 0;
            double bgTotal = 0;
            double elevTotal = 0;
            double bgSampleTotal = sumVector(mBackgroundCounts.getCol(sampleId));

            for (int bucketId = 0; bucketId < bucketCount; ++bucketId)
            {
                double sbCount = scData[bucketId][sampleId];

                boolean isAllocated = false;
                for(Integer bgIndex: bucketGroups)
                {
                    BucketGroup bucketGroup = mBucketGroups.get(bgIndex);

                    if(bucketGroup.hasBucket(bucketId))
                    {
                        isAllocated = true;
                        break;
                    }
                }

                boolean isElevated = elevatedBuckets != null && elevatedBuckets.contains(bucketId);

                if(!isElevated && !isAllocated)
                {
                    bgData[bucketId][sampleId] = sbCount;
                    bgTotal += sbCount;
                    continue;
                }

                if(isElevated)
                    ++lowProbCount;

                ++elevatedCount;

                // split the count between them
                double bgCount = bgData[bucketId][sampleId];

                if(sbCount >= bgCount)
                {
                    double elevCount = sbCount - bgCount;
                    evData[bucketId][sampleId] = elevCount;
                    elevTotal += elevCount;
                    bgTotal += bgCount;
                }
                else
                {
                    // allocate all to elevated- should be rare??
                    bgData[bucketId][sampleId] = 0;
                    evData[bucketId][sampleId] = sbCount;
                    ++expAboveActual;
                    expAboveActualTotal += sbCount;
                    elevTotal += sbCount;

                    if(isElevated)
                        ++lowProbExpAboveActual;
                }
            }

            if(expAboveActual > 0.5 * elevatedCount || expAboveActualTotal > 0.1 * bgSampleTotal || lowProbExpAboveActual > 0)
            {
                LOGGER.debug(String.format("sample(%d) totals(bg=%.0f elev=%.0f bgExp=%.0f act=%.0f) expAboveActual(%d total=%.0f) elevated(%d lowProb=%d clash=%d)",
                        sampleId, bgTotal, elevTotal, bgSampleTotal, mSampleTotals[sampleId],
                        expAboveActual, expAboveActualTotal, elevatedCount, lowProbCount, lowProbExpAboveActual));
            }
        }

        writeSampleMatrixData(elevatedData, "ba_elevated_sc.csv");
        writeSampleMatrixData(mBackgroundCounts, "ba_background_sc.csv");
    }

    private void createSignatures()
    {
        int bucketCount = mSampleCounts.Rows;

        int proposedSigCount = min(mMaxProposedSigs, mBucketGroups.size());

        mProposedSigs = new NmfMatrix(bucketCount, proposedSigCount);
        double[][] sigData = mProposedSigs.getData();
        final double[][] scData = mSampleCounts.getData();

        for(int sigId = 0; sigId < proposedSigCount; ++sigId)
        {
            final BucketGroup bucketGroup = mBucketGroups.get(sigId);

            final List<Integer> sampleIds = bucketGroup.getSampleIds();
            final List<Integer> bucketIds = bucketGroup.getBucketIds();

            // recalc the bucket count totals since merging can double-count some samples
            double[] bucketCounts = new double[bucketCount];
            for(Integer bucketId : bucketIds)
            {
                for(Integer sampleId : sampleIds)
                {
                    bucketCounts[bucketId] += scData[bucketId][sampleId];
                }
            }

            double totalCount = sumVector(bucketCounts);

            // now normalise the signature to percents
            for(int i = 0; i < bucketCount; ++i)
            {
                sigData[i][sigId] = bucketCounts[i] / totalCount;
            }
        }

        writeSignatures(mProposedSigs);
    }

    private boolean areAnySampleBucketsAllocated(int sample, final List<Integer> bucketSubset)
    {
        for (BucketGroup bucketGroup : mBucketGroups)
        {
            if(!bucketGroup.hasSample(sample))
                continue;

            for(Integer bucket : bucketSubset)
            {
                if(bucketGroup.hasBucket(bucket))
                    return true;
            }
        }

        return false;
    }

    private List<Integer> getSampleBucketGroups(int sample)
    {
        List<Integer> groupIds = Lists.newArrayList();

        for (int i = 0; i < mBucketGroups.size(); ++i)
        {
            BucketGroup bucketGroup = mBucketGroups.get(i);

            if (bucketGroup.hasSample(sample))
                groupIds.add(i);
        }

        return groupIds;
    }

    private int getSampleGroupMembershipCount(int sample)
    {
        return getSampleBucketGroups(sample).size();
    }

    private boolean areMostSampleBucketsAllocated(int sample, final List<Integer> sampleBuckets, double reqPercent)
    {
        int allocatedCount = getSampleBucketsAllocated(sample, sampleBuckets);
        return allocatedCount / (double)sampleBuckets.size() >= reqPercent;
    }

    private int getSampleBucketsAllocated(int sample, final List<Integer> sampleBuckets)
    {
        List<Integer> allocatedBuckets = Lists.newArrayList();
        for (BucketGroup bucketGroup : mBucketGroups)
        {
            if (!bucketGroup.hasSample(sample))
                continue;

            for(Integer samBucket : sampleBuckets)
            {
                if(allocatedBuckets.contains(samBucket))
                    continue;

                if (bucketGroup.hasBucket(samBucket))
                    allocatedBuckets.add(samBucket);

                if(allocatedBuckets.size() == sampleBuckets.size())
                    break;
            }

            if(allocatedBuckets.size() == sampleBuckets.size())
                break;
        }

        return allocatedBuckets.size();
    }

    private double[] extractBucketCountSubset(int sam1, final List<Integer> bucketSubset)
    {
        // extract the counts for the specified subset, leaving the rest zeroed
        double[] vec = new double[mSampleCounts.Rows];

        final double[][] scData = mSampleCounts.getData();

        for(int i = 0; i < bucketSubset.size(); ++i)
        {
            int bucketIndex = bucketSubset.get(i);
            vec[bucketIndex] = scData[bucketIndex][sam1];
        }

        return vec;
    }

    private double calcSharedCSS(final double[] set1, final double[] set2)
    {
        return calcCSS(set1, set2, true);
    }

    private final String getSampleName(int sampleId)
    {
        return mDataCollection.getFieldNames().get(sampleId);
    }

    private HashMap<String,Integer> populateSampleCategoryMap(List<Integer> sampleIds)
    {
        final List<String> fieldNames = mExtSampleData.getFieldNames();
        final List<String> sampleNames = mDataCollection.getFieldNames();

        int categoryCount = fieldNames.size() - 1; // since sampleId is the first column

        HashMap<String,Integer> categoryCounts = new HashMap();

        for(Integer sampleId : sampleIds)
        {
            final String sampleName = sampleNames.get(sampleId);

            // now locate all info about this sample from the file
            List<String> sampleData = getSampleExtData(sampleName);

            if(sampleData == null)
                continue;

            for(int c = 1; c <= categoryCount; ++c) // since field name is first
            {
                final String catValue = sampleData.get(c);

                if (catValue.isEmpty())
                    continue;

                String catName = fieldNames.get(c);

                if (c <= CATEGORY_COL_COUNT)
                    catName += CATEGORY_DELIM + catValue;

                if (!categoryCounts.containsKey(catName))
                    categoryCounts.put(catName, 1);
                else
                    categoryCounts.put(catName, categoryCounts.get(catName) + 1);
            }
        }

        return categoryCounts;
    }

    private String CATEGORY_DELIM = "-";

    private String extractCategoryName(String categoryVal)
    {
        if(!categoryVal.contains(CATEGORY_DELIM))
            return "";

        return categoryVal.substring(0, categoryVal.indexOf(CATEGORY_DELIM));
    }

    private String extractCategoryValue(String categoryVal)
    {
        if(!categoryVal.contains(CATEGORY_DELIM))
            return categoryVal;

        return categoryVal.substring(categoryVal.indexOf(CATEGORY_DELIM)+1);
    }

    private final List<String> getSampleExtData(final String sampleName)
    {
        final List<List<String>> extSampleData = mExtSampleData.getStringData();

        for(final List<String> sampleData : extSampleData)
        {
            if(sampleData.get(COL_SAMPLE_ID).equals(sampleName))
                return sampleData;
        }

        return null;
    }

    private void writeSampleData()
    {
        try
        {
            BufferedWriter writer = getNewFile(mOutputDir,mOutputFileId + "_ba_sample_alloc.csv");

            writer.write("SampleIndex,SampleId,Status,BgIds,BgCancerType,BgEffects,Buckets");
            writer.newLine();

            for(int i = 0; i < mSampleCounts.Cols; ++i)
            {
                final List<Integer> sampleBuckets = mAllSampleBucketGroups.get(i);

                if(sampleBuckets == null)
                    continue;;

                final List<Integer> bgIds = getSampleBucketGroups(i);
                final String sampleName = getSampleName(i);

                int bucketsAllocated = getSampleBucketsAllocated(i, sampleBuckets);

                String status = "Unalloc";
                if(bucketsAllocated/(double)sampleBuckets.size() >= 0.9)
                {
                    status = "Alloc";
                }
                else if(bucketsAllocated > 0)
                {
                    status = "Partial";
                }

                writer.write(String.format("%d,%s,%s", i, sampleName, status));

                if(!bgIds.isEmpty())
                {
                    final BucketGroup bucketGroup = mBucketGroups.get(bgIds.get(0));
                    String bgIdsStr = bgIds.toString().substring(1, bgIds.toString().length()-1);
                    if(bgIds.size() > 1)
                        bgIdsStr = bgIdsStr.replaceAll(", ", ";");

                    writer.write(String.format(",%s,%s,%s",
                            bgIdsStr, bucketGroup.getCancerType(), bucketGroup.getEffects()));
                }
                else
                {
                    writer.write(",,,");
                }

                String bucketIdsStr = sampleBuckets.toString().substring(1, sampleBuckets.toString().length()-1); // remove []
                if(sampleBuckets.size() > 1)
                    bucketIdsStr = bucketIdsStr.replaceAll(", ", ";");

                writer.write(String.format(",%s", bucketIdsStr));
                writer.newLine();
            }

            writer.close();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile");
        }
    }

    public void writeSignatures(final NmfMatrix signatures)
    {
        try
        {
            BufferedWriter writer = getNewFile(mOutputDir,mOutputFileId + "_ba_sigs.csv");

            int i = 0;
            for(; i < signatures.Cols-1; ++i)
            {
                writer.write(String.format("%d,", i));
            }
            writer.write(String.format("%d", i));

            writer.newLine();

            writeMatrixData(writer, signatures, false);

            writer.close();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile");
        }
    }

    public void writeSampleMatrixData(final NmfMatrix sampleData, final String fileId)
    {
        try
        {
            BufferedWriter writer = getNewFile(mOutputDir, mOutputFileId + fileId); // eg "_ba_contribs.csv"

            List<String> sampleNames = mDataCollection.getFieldNames();

            int i = 0;
            for(; i < sampleNames.size()-1; ++i)
            {
                writer.write(String.format("%s,", sampleNames.get(i)));
            }
            writer.write(String.format("%s", sampleNames.get(i)));

            writer.newLine();

            writeMatrixData(writer, sampleData, false);

            writer.close();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile");
        }
    }

}
