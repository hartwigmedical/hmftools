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
import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.getTopCssPairs;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getNewFile;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sumVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.writeMatrixData;
import static com.hartwig.hmftools.data_analyser.calcs.NmfConfig.NMF_FS_MIN_SAMPLES;
import static com.hartwig.hmftools.data_analyser.calcs.NmfConfig.NMF_REF_SIG_FILE;
import static com.hartwig.hmftools.data_analyser.types.BucketGroup.getCombinedBuckets;
import static com.hartwig.hmftools.data_analyser.types.BucketGroup.getMatchingBucketList;
import static com.hartwig.hmftools.data_analyser.types.BucketPair.RATIO_MAX;
import static com.hartwig.hmftools.data_analyser.types.BucketPair.RATIO_MIN;
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
import com.hartwig.hmftools.data_analyser.types.BucketPair;
import com.hartwig.hmftools.data_analyser.types.GenericDataCollection;
import com.hartwig.hmftools.data_analyser.types.NmfMatrix;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class BucketAnalyser {

    private static final Logger LOGGER = LogManager.getLogger(NmfManager.class);

    private String mOutputDir;
    private String mOutputFileId;

    private GenericDataCollection mDataCollection;

    private NmfMatrix mSampleCounts;

    private NmfMatrix mBucketMedianRatios;
    private NmfMatrix mBucketProbs;
    private double[] mSampleTotals;
    HashMap<Integer, List<Integer>> mSampleBucketGroups;
    private NmfMatrix mReferenceSigs;

    // external sample data for verification and correlation
    GenericDataCollection mExtSampleData;
    HashMap<String,Integer> mExtCategoriesMap;

    private List<BucketGroup> mBucketGroups;

    private static int MIN_BG_BUCKET_COUNT = 20; // below which a bucket ratio is considered invalid
    private static double BUCKET_COUNTS_CSS_THRESHOLD = 0.995;
    private static double BUCKET_MERGE_CSS_THRESHOLD = 0.99; // lower since is operated over potentially non-matching bucket sets
    private static double MIN_BUCKET_ALLOCATION = 0.9; // consider a sample fully allocated if X% of its buckets are attributed to groups
    private static double DOMINANT_CATEGORY_PERCENT = 0.80; /// mark a group with a category if X% of samples in it have this attribute (eg cancer type, UV)


    // external data file attributes
    private static String BA_EXT_SAMPLE_DATA_FILE = "ba_ext_data_file";
    private static int COL_SAMPLE_ID = 0;
    private static int COL_CANCER_TYPE = 1;
    private static int CATEGORY_COL_COUNT = 3;
    private static String CATEGORY_CANCER_TYPE = "Cancer";
    private static String EFFECT_TRUE = "TRUE";

    public BucketAnalyser()
    {
        mOutputDir = "";
        mOutputFileId = "";

        mDataCollection = null;
        mSampleCounts = null;
        mBucketMedianRatios = null;
        mSampleTotals = null;
        mExtSampleData = null;
        mExtCategoriesMap = null;
        mReferenceSigs = null;

        mBucketGroups = Lists.newArrayList();
    }

    public static void addCmdLineArgs(Options options) {

        options.addOption(BA_EXT_SAMPLE_DATA_FILE, true, "Sample external data");
    }

    public void initialise(GenericDataCollection collection, final CommandLine cmd)
    {
        mDataCollection = collection;
        mOutputFileId = cmd.getOptionValue(OUTPUT_FILE_ID);
        mOutputDir = cmd.getOptionValue(OUTPUT_DIR);

        double minSamplePerc = cmd.hasOption(NMF_FS_MIN_SAMPLES) ? Double.parseDouble(cmd.getOptionValue(NMF_FS_MIN_SAMPLES)) : 0.01;

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
        formBucketGroups(0.95);

        // now consider samples belonging to more than 1 bucket group
        multiAssignSamples(0.3);

        mergeBucketGroups();

        // try again - is anything more added?
        multiAssignSamples(0.3);

        perfCounter.stop();

        logBucketGroups(false);

        perfCounter.start("AnalyseResults");
        analyseGroupsVsExtData();
        logUnallocatedSamples();
        compareRefSigs();
        perfCounter.stop();

        /*
        // form signatures and contributions from these groups
        perfCounter.start("CreateSigs");
        createSignatures();
        perfCounter.stop();
        */

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

        double bucketMedRange = 0.02;
        int minMedIndex = (int) floor(sampleCount * 0.5 - bucketMedRange);
        int maxMedIndex = (int) ceil(sampleCount * 0.5 + bucketMedRange);
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

    private static double MIN_ELEVATED_RATIO = 2;
    private static double MAX_ELEVATED_PROB = 1e-12;

    private void collectElevatedSampleBuckets()
    {
        int bucketCount = mSampleCounts.Rows;
        int sampleCount = mSampleCounts.Cols;

        mSampleBucketGroups = new HashMap();

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

            mSampleBucketGroups.put(i, bucketList);
            // LOGGER.debug("sample({}) has {} elevated buckets", i, bucketList.size());
        }

        LOGGER.debug(String.format("samples with elevated buckets: count(%d perc=%.2f), buckets(%d perc=%.3f)",
                mSampleBucketGroups.size(), mSampleBucketGroups.size()/(double)sampleCount,
                totalCount, totalCount/(double)(bucketCount*sampleCount)));
    }

    private boolean formBucketGroups(double reqMatchPercent)
    {
        // forms groups out of samples with similarly elevated buckets
        // returns true if groups were created or added to
        int sampleCount = mSampleCounts.Cols;

        int groupsAdjusted = 0;
        int groupsCreated = 0;

        for (int samIndex1 = 0; samIndex1 < sampleCount; ++samIndex1)
        {
            final List<Integer> bl1 = mSampleBucketGroups.get(samIndex1);

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

                    if(bcCss >= BUCKET_COUNTS_CSS_THRESHOLD)
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
                mSampleBucketGroups.remove(samIndex1);
                continue;
            }

            for (int samIndex2 = samIndex1 + 1; samIndex2 < sampleCount; ++samIndex2)
            {
                final List<Integer> bl2 = mSampleBucketGroups.get(samIndex2);

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

                    if(bcCss >= BUCKET_COUNTS_CSS_THRESHOLD)
                    {
                        BucketGroup bucketGroup = new BucketGroup(mBucketGroups.size());
                        bucketGroup.addBuckets(commonBuckets);
                        bucketGroup.addSample(samIndex1, sc1);
                        bucketGroup.addSample(samIndex2, sc2);

                        LOGGER.debug(String.format("added bg(%d) samples(%d and %d) with buckets(s1=%d and s2=%d matched=%d) css(%.4f)",
                                bucketGroup.getId(), samIndex1, samIndex2, bl1.size(), bl2.size(), commonBuckets.size(), bcCss));

                        mBucketGroups.add(bucketGroup);
                        ++groupsCreated;

                        mSampleBucketGroups.remove(samIndex2);
                        mSampleBucketGroups.remove(samIndex1);
                        break;
                    }
                }
            }
        }

        if(mBucketGroups.isEmpty())
        {
            LOGGER.debug("no bucket groups created");
            return false;
        }

        LOGGER.debug("bucket groups created({}) additions({}) total({}), unallocated elevated samples({})",
                groupsCreated, groupsAdjusted, mBucketGroups.size(), mSampleBucketGroups.size());

        // log from most comprehensive to least
        Collections.sort(mBucketGroups);

        return true;
    }

    private boolean multiAssignSamples(double reqMatchPercent)
    {
        LOGGER.debug(String.format("checking multi-assignments for elevated samples(%d) and %d bucket groups",
                mSampleBucketGroups.size(), mBucketGroups.size()));

        // check unallocated samples against each other and existing groups, looking for partial overlap
        int sampleCount = mSampleCounts.Cols;

        int groupsAdjusted = 0;
        int groupsCreated = 0;

        for (int samIndex1 = 0; samIndex1 < sampleCount; ++samIndex1)
        {
            final List<Integer> bl1 = mSampleBucketGroups.get(samIndex1);

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

                final List<Integer> commonBuckets = getMatchingBucketList(bucketGroup.getBucketIds(), bl1);

                if(areSampleBucketsAllocated(samIndex1, commonBuckets))
                {
                    // LOGGER.debug("bg({}) vs sample({}) when already allocated to another BG", bucketGroup.getId(), samIndex1);
                    continue;
                }

                // look for a sample with an excess of elevated buckets
                if (bl1.size() <= commonBuckets.size())
                    continue;

                // test similarity for the common buckets
                double[] sc1 = extractBucketCountSubset(samIndex1, commonBuckets);
                double bcCss = calcSharedCSS(sc1, bucketGroup.getBucketCounts());

                double sampleMatch = commonBuckets.size() / (double)bl1.size();

                if (sampleMatch >= reqMatchPercent && bcCss >= BUCKET_COUNTS_CSS_THRESHOLD)
                {
                    assignmentScores[i] = sampleMatch * bcCss;
                    hasPossibleAssignments = true;

//                    LOGGER.debug(String.format("sample(%d) could match bg(%d) with buckets(grp=%d sam=%d match=%d) totalSamples(%d) css(%.4f)",
//                            samIndex1, bucketGroup.getId(), bucketGroup.getBucketIds().size(), bl1.size(),
//                            commonBuckets.size(), bucketGroup.getSampleIds().size(), bcCss));

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

                    if(areSampleBucketsAllocated(samIndex1, commonBuckets))
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
                    mSampleBucketGroups.remove(samIndex1);
                    continue;
                }
            }

            // otherwise keep checking amongst unallocated samples in a similar manner
            boolean sam1Allocated = false;

            for (int samIndex2 = samIndex1 + 1; samIndex2 < sampleCount; ++samIndex2)
            {
                final List<Integer> bl2 = mSampleBucketGroups.get(samIndex2);

                if (bl2 == null)
                    continue;

                final List<Integer> commonBuckets = getMatchingBucketList(bl1, bl2);

                if(areSampleBucketsAllocated(samIndex1, commonBuckets))
                {
                    continue;
                }

                double[] sc1 = extractBucketCountSubset(samIndex1, commonBuckets);
                double[] sc2 = extractBucketCountSubset(samIndex2, commonBuckets);
                double bcCss = calcSharedCSS(sc1, sc2);

                double sample1Match = commonBuckets.size() / (double)bl1.size();
                double sample2Match = commonBuckets.size() / (double)bl2.size();
                double minMatch = min(sample1Match, sample2Match);

                if(minMatch > 1 && minMatch >= reqMatchPercent && bcCss >= BUCKET_COUNTS_CSS_THRESHOLD)
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
                        mSampleBucketGroups.remove(samIndex2);
                        continue;
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
                mSampleBucketGroups.remove(samIndex1);
            }
            else
            {
                ++samIndex1;
            }
        }

        if(groupsAdjusted == 0 && groupsCreated == 0)
        {
            LOGGER.debug("no new/adjusted bucket groups");
            return false;
        }

        LOGGER.debug("bucket groups created({}) adjusted({}) total({}), unallocated elevated samples({})",
                groupsCreated, groupsAdjusted, mBucketGroups.size(), mSampleBucketGroups.size());

        // log from most comprehensive to least
        Collections.sort(mBucketGroups);

        return true;
    }

    private boolean areSampleBucketsAllocated(int sample, final List<Integer> bucketSubset)
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

    private int getSampleGroupMembershipCount(int sample)
    {
        int groupCount = 0;

        for (BucketGroup bucketGroup : mBucketGroups)
        {
            if (bucketGroup.hasSample(sample))
                ++groupCount;
        }

        return groupCount;
    }

    private boolean areMostSampleBucketsAllocated(int sample, final List<Integer> sampleBuckets, double reqPercent)
    {
        int allocatedCount = getSampleBucketsAllocated(sample, sampleBuckets);
        return allocatedCount / (double)sampleBuckets.size() >= reqPercent;
    }

    private int getSampleBucketsAllocated(int sample, final List<Integer> sampleBuckets)
    {
        int allocatedCount = 0;
        for(Integer samBucket : sampleBuckets)
        {
            boolean isAllocated = false;

            for (BucketGroup bucketGroup : mBucketGroups)
            {
                if (!bucketGroup.hasSample(sample))
                    continue;

                if (bucketGroup.hasBucket(samBucket))
                {
                    isAllocated = true;
                    break;
                }
            }

            if(isAllocated)
                ++allocatedCount;
        }

        return allocatedCount;
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

    private void mergeBucketGroups()
    {
        double reqMatchPercent = 0.85;

        int bgIndex1 = 0;
        while (bgIndex1 < mBucketGroups.size())
        {
            BucketGroup bg1 = mBucketGroups.get(bgIndex1);

            int bgIndex2 = bgIndex1+1;
            while (bgIndex2 < mBucketGroups.size())
            {
                BucketGroup bg2 = mBucketGroups.get(bgIndex2);

                final List<Integer> commonBuckets = getMatchingBucketList(bg1.getBucketIds(), bg2.getBucketIds());

                int commmonBucketCount = commonBuckets.size();
                double bg1Match = commmonBucketCount / (double) bg1.getBucketIds().size();
                double bg2Match = commmonBucketCount / (double) bg2.getBucketIds().size();
                double minMatch = min(bg1Match, bg2Match);
                double bcCss = calcSharedCSS(bg1.getBucketCounts(), bg2.getBucketCounts());

                if (minMatch >= reqMatchPercent && bcCss >= BUCKET_MERGE_CSS_THRESHOLD)
                {
                    int bg1SC = bg1.getSampleIds().size();
                    int bg2SC = bg2.getSampleIds().size();

                    if(minMatch < 1)
                    {
                        // should new buckets be merged in? or just go with the common set or most prevalent set
                        if (bg1SC > 1.5 * bg2SC)
                        {
                            bg1.merge(bg2.getSampleIds(), bg2.getBucketCounts());
                        }
                        else if (bg2SC > 2 * bg1SC)
                        {
                            bg2.merge(bg1.getSampleIds(), bg1.getBucketCounts());
                        }
                        else
                        {
                            // go with the common set only
                            bg1.merge(bg2.getSampleIds(), bg2.getBucketCounts());
                            bg1.reduceToBucketSet(commonBuckets);
                        }
                    }

                    LOGGER.debug(String.format("bg(%d) merge with bg(%d) buckets(bg1=%d bg2=%d match=%d) css(%.4f) samples(bg1=%d bg2=%d total=%d)",
                            bg1.getId(), bg2.getId(), bg1.getBucketIds().size(), bg2.getBucketIds().size(), commmonBucketCount,
                            bcCss, bg1SC, bg2SC, bg1SC + bg2SC));

                    mBucketGroups.remove(bgIndex2);
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

                ++bgIndex2;
            }

            ++bgIndex1;
        }
    }

    private double calcSharedCSS(final double[] set1, final double[] set2)
    {
        return calcCSS(set1, set2, true);
    }

    private void logBucketGroups(boolean clearTypesOnly)
    {
        // log top groups
        int maxToLog = clearTypesOnly ? mBucketGroups.size() : min(mBucketGroups.size(), 40);
        int minScore = clearTypesOnly ? 0 : 10; // 4 samples x 5 buckets, or 2 x 10

        for (int i = 0; i < maxToLog; ++i)
        {
            BucketGroup bucketGroup = mBucketGroups.get(i);

            if (bucketGroup.calcScore() < minScore)
                break;

            if(clearTypesOnly)
            {
                if(bucketGroup.getCancerType().isEmpty())
                    continue;

                LOGGER.debug(String.format("bg(%d) cancer(%s) score(%.0f) samples(%d) buckets(%d: %s) effects(%s)",
                        bucketGroup.getId(), bucketGroup.getCancerType(), bucketGroup.calcScore(), bucketGroup.getSampleIds().size(),
                        bucketGroup.getBucketIds().size(), bucketGroup.getBucketIds().toString(), bucketGroup.getEffects()));
            }
            else
            {
                LOGGER.debug(String.format("bg(%d) score(%.0f) samples(%d) buckets(%d)",
                        bucketGroup.getId(), bucketGroup.calcScore(), bucketGroup.getSampleIds().size(),
                        bucketGroup.getBucketIds().size()));
            }
        }
    }

    private void logUnallocatedSamples()
    {
        final List<String> categories = mExtSampleData.getFieldNames();

        int partiallyGrouped = 0;
        int noMatchCount = 0;
        int someMatchCount = 0;

        for(Map.Entry<Integer, List<Integer>> entry : mSampleBucketGroups.entrySet())
        {
            int sampleId = entry.getKey();
            final List<Integer> samBucketList = entry.getValue();

            final String sampleName = getSampleName(sampleId);

            final List<String> sampleData = getSampleExtData(sampleName);
            if(sampleData == null)
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

            int groupCount = getSampleGroupMembershipCount(sampleId);
            if(groupCount > 0)
            {
                ++partiallyGrouped;
                int bucketsAllocated = getSampleBucketsAllocated(sampleId, samBucketList);
                LOGGER.debug("sample({}: {}) partially allocated to {} groups: buckets({} alloc={}) cancer({}) effects({}: {})",
                        sampleId, sampleName, groupCount, samBucketList.size(), bucketsAllocated, cancerType, effectsCount, effects);
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

        LOGGER.debug("unallocated samples({}) partiallyGrouped({}) partialMatch({}) noMatch({})",
                mSampleBucketGroups.size(), partiallyGrouped, someMatchCount, noMatchCount);
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

        logBucketGroups(true);

        // need to figure out why there are multiple bucket groups for the same underlying condition (be it cancer type or process eg AID/APOBEC)
        // to do this, for each category, log all related buckets and any sample overlapping a decent proportion of the most common buckets
        // for(Map.Entry<String,Integer>)
    }

    private void compareRefSigs()
    {
        if(mReferenceSigs == null)
            return;

        // form a temporary matrix of signature data from the bucket groups, then copmare with extenal sigs
        NmfMatrix possibleSigs = new NmfMatrix(mSampleCounts.Rows, mBucketGroups.size());

        for (int i = 0; i < mBucketGroups.size(); ++i)
        {
            BucketGroup bucketGroup = mBucketGroups.get(i);
            possibleSigs.setCol(i, bucketGroup.getBucketRatios());
        }

        List<double[]> cssResults = getTopCssPairs(mReferenceSigs, possibleSigs, 0.95, true, false);

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

                LOGGER.debug(String.format("external ref sig(%d) matches bg(%d: ct=%s eff=%s) with css(%.4f)",
                        externalSigId, bucketGroup.getId(), bucketGroup.getCancerType(), bucketGroup.getEffects(), result[CSSR_VAL]));
            }
        }
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

    private void createSignatures()
    {
        // convert the highest scoring bucket groups into signatures

        int bucketCount = mSampleCounts.Rows;
        int sampleCount = mSampleCounts.Cols;

        List<Integer> uniqueSampleIds = Lists.newArrayList();

        int maxProposedSigs = 20;
        int proposedSigCount = min(maxProposedSigs, mBucketGroups.size());

        NmfMatrix proposedSigs = new NmfMatrix(bucketCount, proposedSigCount);
        double[][] sigData = proposedSigs.getData();
        final double[][] scData = mSampleCounts.getData();

        NmfMatrix contributions = new NmfMatrix(proposedSigCount, sampleCount); // uniqueSampleIds.size()
        double[][] contribData = contributions.getData();

        for(int sigId = 0; sigId < proposedSigCount; ++sigId)
        {
            final BucketGroup bucketGroup = mBucketGroups.get(sigId);

            final List<Integer> sampleIds = bucketGroup.getSampleIds();

            double totalCount = 0;

            for(Integer sampleId : sampleIds)
            {
                uniqueSampleIds.add(sampleId);

                for(int j = 0; j < bucketCount; ++j)
                {
                    if(!bucketGroup.hasBucket(j))
                        continue;

                    double sbCount = scData[j][sampleId];
                    sigData[j][sigId] += sbCount;
                    totalCount += sbCount;

                    // each sample will have their buckets from this group allocated to this signature
                    contribData[sigId][sampleId] += sbCount;
                }
            }

            // now normalise the signature to percents
            for(int j = 0; j < bucketCount; ++j)
            {
                sigData[j][sigId] = sigData[j][sigId] / totalCount;
            }
        }

//        for(int sigId = 0; sigId < proposedSigCount; ++sigId)
//        {
//            final BucketGroup bucketGroup = mBucketGroups.get(sigId);
//
//            for (int i = 0; i < sampleCount; ++i)
//            {
//                if (bucketGroup.hasSample(i))
//                    contribData[sigId][i] = 1;
//            }
//        }

        writeSignatures(proposedSigs);
        writeContributions(contributions);
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

    public void writeContributions(final NmfMatrix contributions)
    {
        try
        {
            BufferedWriter writer = getNewFile(mOutputDir, mOutputFileId + "_ba_contribs.csv");

            List<String> sampleNames = mDataCollection.getFieldNames();

            int i = 0;
            for(; i < sampleNames.size()-1; ++i)
            {
                writer.write(String.format("%s,", sampleNames.get(i)));
            }
            writer.write(String.format("%s", sampleNames.get(i)));

            writer.newLine();

            writeMatrixData(writer, contributions, false);

            writer.close();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile");
        }
    }

}
