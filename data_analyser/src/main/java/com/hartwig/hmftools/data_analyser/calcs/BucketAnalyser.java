package com.hartwig.hmftools.data_analyser.calcs;

import static java.lang.Math.abs;
import static java.lang.Math.log;
import static java.lang.Math.log10;
import static java.lang.Math.max;
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
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.calcBestFitWithinProbability;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.calcLinearLeastSquares;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getCombinedList;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getDiffList;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getMatchingList;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getNewFile;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.listToArray;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sizeToStr;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sumVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.writeMatrixData;
import static com.hartwig.hmftools.data_analyser.calcs.NmfConfig.NMF_REF_SIG_FILE;
import static com.hartwig.hmftools.data_analyser.types.GenericDataCollection.GD_TYPE_STRING;
import static com.hartwig.hmftools.data_analyser.types.NmfMatrix.redimension;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.PerformanceCounter;
import com.hartwig.hmftools.data_analyser.loaders.GenericDataLoader;
import com.hartwig.hmftools.data_analyser.types.BucketFamily;
import com.hartwig.hmftools.data_analyser.types.BucketGroup;
import com.hartwig.hmftools.data_analyser.types.GenericDataCollection;
import com.hartwig.hmftools.data_analyser.types.NmfMatrix;
import com.hartwig.hmftools.data_analyser.types.SampleData;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class BucketAnalyser {

    private static final Logger LOGGER = LogManager.getLogger(NmfManager.class);

    private GenericDataCollection mDataCollection;

    private NmfMatrix mSampleCounts;

    // for convenience
    private double[] mSampleTotals;
    private double mTotalCount;
    private int mBucketCount;
    private int mSampleCount;

    private Map<String, List<Double>> mBucketMediansMap; // cancer-type to median bucket ratios (ie background sigs)
    private NmfMatrix mBucketMedianRatios;
    private NmfMatrix mBucketProbs;
    private NmfMatrix mBackgroundCounts;
    private NmfMatrix mElevatedCounts; // actual - expected, capped at zero
    private int mElevatedCount;
    private NmfMatrix mPermittedElevRange;
    private NmfMatrix mPermittedBgRange;
    private List<Double> mSampleBgAllocations;

    HashMap<Integer, List<Integer>> mAllSampleBucketGroups;
    private NmfMatrix mProposedSigs;
    private List<Integer> mSigToBgMapping;
    private NmfMatrix mReferenceSigs;

    // external sample data for verification and correlation
    private List<SampleData> mSampleData;
    private GenericDataCollection mExtSampleData;
    private HashMap<String,Integer> mExtCategoriesMap;
    private HashMap<String, List<Integer>> mCancerSamplesMap;
    private boolean mLoadedSampleCalcData;

    private List<BucketGroup> mBucketGroups;
    private List<BucketGroup> mBackgroundGroups;
    private List<BucketFamily> mBucketFamilies;
    private boolean mBucketGroupsLinked;

    // config
    private String mOutputDir;
    private String mOutputFileId;
    private double mHighCssThreshold; // CSS level for samples or groups to be consider similar
    private double mHighRequiredMatch; // high-level required number of buckets to match between samples or groups
    private int mMaxProposedSigs;
    private double mSigCssThreshold; // avoid creating sigs with similarity lower than this
    private int mMutationalLoadCap;

    private static String BA_CSS_HIGH_THRESHOLD = "ba_css_high";
    private static String BA_MAX_PROPOSED_SIGS = "ba_max_proposed_sigs";
    private static String BA_CSS_SIG_THRESHOLD = "ba_css_proposed_sigs";

    // constraint constants - consider moing any allocation-related ones to config to make them visible
    private static double MIN_BUCKET_ALLOCATION = 0.9; // consider a sample fully allocated if X% of its buckets are attributed to groups
    private static double DOMINANT_CATEGORY_PERCENT = 0.7; /// mark a group with a category if X% of samples in it have this attribute (eg cancer type, UV)
    private static double MIN_ELEVATION_RATIO = 1.2;
    private static double MAX_ELEVATED_PROB = 1e-12;
    private int MIN_BUCKET_COUNT_OVERLAP = 5; // used for pairings of samples with reduced bucket overlap
    private static double PERMITTED_PROB_NOISE = 1e-4;

    // external data file attributes
    private static String BA_EXT_SAMPLE_DATA_FILE = "ba_ext_data_file";
    private static int COL_SAMPLE_ID = 0;
    private static int COL_CANCER_TYPE = 1;
    private static int CATEGORY_COL_COUNT = 3;
    private static String CATEGORY_CANCER_TYPE = "Cancer";
    private static String EFFECT_TRUE = "TRUE";
    private static String BA_SAMPLE_CALC_DATA_FILE = "ba_sam_calc_data_file";

    private static int SCD_COL_SAMPLE_ID = 0;
    private static int SCD_COL_SAMPLE_NAME = 1;
    private static int SCD_COL_BG_ALLOC = 2;

    private List<Integer> mSampleWatchList;

    public BucketAnalyser()
    {
        mOutputDir = "";
        mOutputFileId = "";

        mDataCollection = null;
        mSampleCounts = null;
        mBucketMedianRatios = null;
        mSampleTotals = null;
        mTotalCount = 0;
        mBucketCount = 0;
        mSampleCount = 0;
        mBackgroundCounts = null;
        mElevatedCounts = null;
        mElevatedCount = 0;
        mPermittedElevRange = null;
        mPermittedBgRange = null;
        mBucketMediansMap = null;
        mSampleBgAllocations = Lists.newArrayList();

        mSampleData = Lists.newArrayList();

        mExtSampleData = null;
        mExtCategoriesMap = null;
        mCancerSamplesMap = null;
        mReferenceSigs = null;
        mProposedSigs = null;
        mSigToBgMapping = Lists.newArrayList();
        mLoadedSampleCalcData = false;

        mBucketGroups = Lists.newArrayList();
        mBackgroundGroups = Lists.newArrayList();
        mBucketFamilies = Lists.newArrayList();
        mBucketGroupsLinked = false;

        mHighCssThreshold = 0.995;
        mHighRequiredMatch = 0.95;
        mMaxProposedSigs = 0;
        mSigCssThreshold = mHighCssThreshold * 0.95;
        mMutationalLoadCap = 0;

        mSampleWatchList = Lists.newArrayList();

//        mSampleWatchList.add(1424);
//        mSampleWatchList.add(2024);
//        mSampleWatchList.add(2617);
    }

    public static void addCmdLineArgs(Options options)
    {

        options.addOption(BA_EXT_SAMPLE_DATA_FILE, true, "Sample external data");
        options.addOption(BA_CSS_HIGH_THRESHOLD, true, "Cosine sim for high-match test");
        options.addOption(BA_CSS_SIG_THRESHOLD, true, "Cosine sim for comparing proposed sigs");
        options.addOption(BA_MAX_PROPOSED_SIGS, true, "Maximum number of bucket groups to turn into proposed sigs");
        options.addOption(BA_SAMPLE_CALC_DATA_FILE, true, "Optional: file containing computed data per sample");
    }

    public boolean initialise(GenericDataCollection collection, final CommandLine cmd)
    {
        mDataCollection = collection;
        mOutputFileId = cmd.getOptionValue(OUTPUT_FILE_ID);
        mOutputDir = cmd.getOptionValue(OUTPUT_DIR);

        mHighCssThreshold = Double.parseDouble(cmd.getOptionValue(BA_CSS_HIGH_THRESHOLD, "0.995"));
        mHighRequiredMatch = 0.95;
        mMaxProposedSigs = Integer.parseInt(cmd.getOptionValue(BA_MAX_PROPOSED_SIGS, "0"));
        mMutationalLoadCap = 10000;

        LOGGER.info("config: cssThreshold({}) reqMatch({})",mHighCssThreshold, mHighRequiredMatch);

        mSampleCounts = DataUtils.createMatrixFromListData(mDataCollection.getData());
        mSampleCounts.cacheTranspose();
        mSampleCount = mSampleCounts.Cols;
        mBucketCount = mSampleCounts.Rows;

        if(cmd.hasOption(NMF_REF_SIG_FILE))
        {
            GenericDataCollection dataCollection = GenericDataLoader.loadFile(cmd.getOptionValue(NMF_REF_SIG_FILE));
            mReferenceSigs = DataUtils.createMatrixFromListData(dataCollection.getData());
            mReferenceSigs.cacheTranspose();
        }

        LOGGER.info("bucketCount({}) sampleCount({})", mBucketCount, mSampleCount);

        if(cmd.hasOption(BA_EXT_SAMPLE_DATA_FILE))
        {
            final String fileName = cmd.getOptionValue(BA_EXT_SAMPLE_DATA_FILE);
            mExtSampleData = GenericDataLoader.loadFile(fileName, GD_TYPE_STRING);

            List<Integer> sampleIds = Lists.newArrayList();
            for (int i = 0; i < mSampleCount; ++i)
            {
                sampleIds.add(i);
            }

            mExtCategoriesMap = populateSampleCategoryMap(sampleIds);

            LOGGER.debug("registered {} cancer types and sub-categories", mExtCategoriesMap.size());
        }
        else
        {
            // could auto-generate a default collection
            LOGGER.error("sample-cancer type file required");
            return false;
        }

        for(int sampleId = 0; sampleId < mSampleCount; ++sampleId)
        {
            SampleData sample = new SampleData(sampleId);
            sample.setBucketCounts(mSampleCounts.getCol(sampleId));

            if(mExtCategoriesMap != null)
            {
                sample.setSampleName(mDataCollection.getFieldNames().get(sampleId));

                final List<String> extSampleData = getSampleExtData(sample.getSampleName());
                if(extSampleData != null)
                {
                    sample.setCancerType(extSampleData.get(COL_CANCER_TYPE));
                    sample.setCategoryData(extSampleData);
                }
            }

            mSampleData.add(sample);
        }

        populateCancerSamplesMap();

        if(cmd.hasOption(BA_SAMPLE_CALC_DATA_FILE))
        {
            final String fileName = cmd.getOptionValue(BA_SAMPLE_CALC_DATA_FILE);
            GenericDataCollection bgSampleAllocations = GenericDataLoader.loadFile(fileName, GD_TYPE_STRING);
            List<List<String>> dataSet = bgSampleAllocations.getStringData();

            for(final List<String> sampleData : dataSet)
            {
                // first 2 fields are sampleId and sampleName, other fields are data related
                if(sampleData.size() < 3)
                    continue;

                int sampleId = Integer.parseInt(sampleData.get(SCD_COL_SAMPLE_ID));
                double bgAllocation = Double.parseDouble(sampleData.get(SCD_COL_BG_ALLOC));
                mSampleBgAllocations.add(sampleId, bgAllocation);
            }

            LOGGER.debug("loaded {} sample calc data fields", dataSet.size());
            mLoadedSampleCalcData = true;
        }

        return true;
    }

    public void run()
    {
        PerformanceCounter perfCounter = new PerformanceCounter("BucketMeanRatios");

        perfCounter.start("SplitCounts");
        calcBucketMedianData();
        calcSampleBackgroundCounts();
        splitSampleCounts();
        collectElevatedSampleBuckets();
        perfCounter.stop();

//        perfCounter.start("ReportSamplePairs");
//        calcSampleOverlaps();
//        perfCounter.stop();

        // back-ground groups using expected counts
        perfCounter.start("BackgroundGroups");
        formBackgroundBucketGroups();
        // logBucketGroups(false, true);
        perfCounter.stop();

        perfCounter.start("BucketGroups");

        // start by looking for near-exact matches in buckets across samples, forming bucket groups
        formExactBucketGroups();
        //analyseGroupsVsExtData(true, false);
        //logBucketGroups(true, true);

        formBucketGroupsFromSamplePairSubsets();

        allocateSamplesToBucketGroups();

        // broadenBucketGroupDefinitions();

        // checkAllSamplesAgainstBucketGroups();

        // formBucketGroupsFromSubsets();

        mergeBucketGroups();
        collapseBucketSamples();
        recalcBucketData();

        perfCounter.stop();

        logBucketGroups(true, false);

        perfCounter.start("AnalyseResults");
        analyseGroupsVsExtData(true, true);
        formBucketFamilies(true);
        logBucketGroups(true, true);
        // writeBucketGroupOverlap();
        logSampleResults();
        perfCounter.stop();

//        createSignatures();
//        compareSignatures();

        writeBucketGroups();
        writeSampleData();
        writeBackgroundSigs();
        writeSampleMatrixData(mElevatedCounts, "_ba_elevated_counts.csv");
        writeSampleCalcData();

        perfCounter.logStats();
    }

    private void calcBucketMedianData()
    {
        mSampleTotals = new double[mSampleCount];

        for (int i = 0; i < mSampleCount; ++i)
        {
            mSampleTotals[i] = sumVector(mSampleCounts.getCol(i));
        }

        mTotalCount = sumVector(mSampleTotals);

        mBucketMediansMap = new HashMap();

        double bucketMedRange = 0.005;

        double[][] scData = mSampleCounts.getData();

        for (Map.Entry<String, List<Integer>> entry : mCancerSamplesMap.entrySet())
        {
            final String cancerType = entry.getKey();
            final List<Integer> sampleIds = entry.getValue();

            // add the counts from every sample with mutational load below the threshold
            // and then infer a background signature (literally bucket ratios) from those only
            int samplesIncluded = 0;

            List<List<Double>> sampleBucketRatios = Lists.newArrayList();

            for (int j = 0; j < sampleIds.size(); ++j)
            {
                int sampleId = sampleIds.get(j);

                double sampleTotal = mSampleTotals[sampleId];

                if (sampleTotal > mMutationalLoadCap)
                    continue;

                ++samplesIncluded;

                List<Double> bucketRatios = Lists.newArrayList();

                for (int i = 0; i < mBucketCount; ++i)
                {
                    bucketRatios.add(scData[i][sampleId] / sampleTotal);
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
                    List<Double> sampleBucketRatio = sampleBucketRatios.get(j);
                    bucketRatios[j] = sampleBucketRatio.get(i);
                }

                // sort these and then take the median
                List<Integer> sortedRatioIndices = getSortedVectorIndices(bucketRatios, true);
                int medianRatioIndex = sortedRatioIndices.get(medianIndex);
                double medianRatio = bucketRatios[medianRatioIndex];
                medianBucketRatios[i] = medianRatio;
            }

            // convert to a percent
            List<Double> medianRatios = Lists.newArrayList();

            double ratioTotal = sumVector(medianBucketRatios);

            for (int i = 0; i < mBucketCount; ++i)
            {
                medianRatios.add(medianBucketRatios[i] / ratioTotal);
            }

            /*

            int minMedIndex = (int) floor(sampleIds.size() * (0.5 - bucketMedRange));
            int maxMedIndex = (int) ceil(sampleIds.size() * (0.5 + bucketMedRange));

            LOGGER.debug("cancerType({}) calculating median counts from indices({} -> {})", cancerType, minMedIndex, maxMedIndex);

            for (int i = 0; i < mBucketCount; ++i)
            {
                double[] bucketCounts = new double[sampleIds.size()];
                for (int j = 0; j < sampleIds.size(); ++j)
                {
                    int sampleId = sampleIds.get(j);
                    bucketCounts[j] = mSampleCounts.get(i, sampleId);
                }

                final List<Integer> bcSorted = DataUtils.getSortedVectorIndices(bucketCounts, true);

                double countTotal = 0;
                for (int j = minMedIndex; j <= maxMedIndex; ++j)
                {
                    countTotal += bucketCounts[bcSorted.get(j)];
                }

                double medBucketCount = countTotal / (maxMedIndex - minMedIndex + 1);

                // LOGGER.debug(String.format("bucket(%d) median count(%.0f)", i, medBucketCount));

                medianCounts.add(medBucketCount);
            }

            */

            mBucketMediansMap.put(cancerType, medianRatios);
        }
    }

    private void calcSampleBackgroundCounts()
    {
        if(!mSampleBgAllocations.isEmpty())
            return;

        LOGGER.debug("calculating sample background counts");

        final double[][] scData = mSampleCounts.getData();

        for (Map.Entry<String, List<Double>> entry : mBucketMediansMap.entrySet())
        {
            final String type = entry.getKey();

            List<Double> medianRatios = entry.getValue();

            List<Integer> sampleIds = mCancerSamplesMap.get(type);

            final double[] bucketRatios = listToArray(medianRatios);

            for (Integer j : sampleIds)
            {
                double[] sampleBgCounts = new double[mBucketCount];

                for (int i = 0; i < mBucketCount; ++i)
                {
                    int expectedCount = (int) round(bucketRatios[i] * min(mSampleTotals[j], mMutationalLoadCap));
                    sampleBgCounts[i] = min(expectedCount, scData[i][j]);
                }

                double optimalBgCount = calcBestFitWithinProbability(bucketRatios, sampleBgCounts, 0.99, 0.01);
                mSampleBgAllocations.add(optimalBgCount);
            }
        }
    }

    private void splitSampleCounts()
    {
        LOGGER.debug("splitting sample counts");

        // work out bucket median values (literally 50th percentile values
        mBucketMedianRatios = new NmfMatrix(mBucketCount, mSampleCount);
        mBucketProbs = new NmfMatrix(mBucketCount, mSampleCount);
        mBackgroundCounts = new NmfMatrix(mBucketCount, mSampleCount);
        mElevatedCounts = new NmfMatrix(mBucketCount, mSampleCount);
        mPermittedElevRange = new NmfMatrix(mBucketCount, mSampleCount);
        mPermittedBgRange = new NmfMatrix(mBucketCount, mSampleCount);

        double[][] brData = mBucketMedianRatios.getData();
        double[][] probData = mBucketProbs.getData();
        double[][] bgData = mBackgroundCounts.getData();
        double[][] elevData = mElevatedCounts.getData();
        double[][] permBgRangeData = mPermittedBgRange.getData();
        double[][] permElevRangeData = mPermittedElevRange.getData();
        final double[][] scData = mSampleCounts.getData();

        Map<Integer, Integer> rangeMap = new HashMap(); // cache range values

        int probExpMax = 15;
        int[] probFrequencies = new int[probExpMax+1];
        double minProb = pow(10, -probExpMax);
        int zeroProbIndex = probFrequencies.length-1;

        int gridSize = mSampleCount * mBucketCount;

        for(Map.Entry<String, List<Double>> entry: mBucketMediansMap.entrySet())
        {
            final String type = entry.getKey();

            List<Double> medianRatios = entry.getValue();

            List<Integer> sampleIds = mCancerSamplesMap.get(type);

            for (int i = 0; i < mBucketCount; ++i)
            {
                double bucketMedianRatio = medianRatios.get(i);

                for (int samIndex = 0; samIndex < sampleIds.size(); ++samIndex)
                {
                    int j = sampleIds.get(samIndex);

                    double bgAllocation = mSampleBgAllocations.get(samIndex);

                    int sbCount = (int)scData[i][j];

                    // int expectedCount = (int) round(bucketMedianRatio * min(mSampleTotals[j], mMutationalLoadCap));
                    int backgroundCount = (int)round(bucketMedianRatio * bgAllocation);

                    // expData[i][j] = expectedCount;
                    // if(expectedCount < backgroundCount)
                    //    backgroundCount = expectedCount;

                    if(sbCount < backgroundCount) // must be capped by the actual count
                        backgroundCount = sbCount;

                    bgData[i][j] = backgroundCount;

                    if (backgroundCount > 0)
                    {
                        Integer rangeVal = rangeMap.get(backgroundCount);
                        if (rangeVal == null)
                        {
                            rangeVal = CosineSim.calcPoissonRangeGivenProb(backgroundCount, PERMITTED_PROB_NOISE);
                            permBgRangeData[i][j] = rangeVal;
                            rangeMap.put(backgroundCount, rangeVal);
                        }
                        else
                        {
                            permBgRangeData[i][j] = rangeVal;
                        }

                        brData[i][j] = sbCount / (double)backgroundCount;
                    }

                    int elevatedCount = sbCount - backgroundCount;

                    // compute a range for Poisson noise around this elevated count
                    if (elevatedCount > 0)
                    {
                        elevData[i][j] = elevatedCount;
                        mElevatedCount += elevatedCount;

                        Integer rangeVal = rangeMap.get(elevatedCount);
                        if (rangeVal == null)
                        {
                            rangeVal = CosineSim.calcPoissonRangeGivenProb(elevatedCount, PERMITTED_PROB_NOISE);
                            permElevRangeData[i][j] = rangeVal;
                            rangeMap.put(elevatedCount, rangeVal);
                        }
                        else
                        {
                            permElevRangeData[i][j] = rangeVal;
                        }
                    }

                    double prob = 1;

                    if (backgroundCount == 0)
                    {
                        prob = 0;
                    }
                    else if (elevatedCount > 0)
                    {
                        PoissonDistribution poisDist = new PoissonDistribution(backgroundCount);
                        prob = 1 - poisDist.cumulativeProbability(sbCount - 1);
                        prob = min(prob * gridSize, 1); // apply false discovery rate being # tests
                    }

                    probData[i][j] = prob;

                    if (prob > minProb && prob < 0.1)
                    {
                        int baseProb = -(int) round(log10(prob));

                        if (baseProb >= 0 && baseProb <= probExpMax)
                            probFrequencies[baseProb] += 1;
                    }
                    else if (prob < minProb)
                    {
                        // allocate to the last slot
                        probFrequencies[zeroProbIndex] += 1;
                    }
                }
            }
        }

        mBackgroundCounts.cacheTranspose();
        mElevatedCounts.cacheTranspose();

        for (int i = 0; i < probFrequencies.length-1; ++i)
        {
            LOGGER.debug(String.format("probability(1e-%d) freq(%d) percOfTotal(%.4f)",
                    i, probFrequencies[i], probFrequencies[i]/(double)gridSize));
        }

        LOGGER.debug(String.format("probability(zero) freq(%d) percOfTotal(%.4f)",
                probFrequencies[zeroProbIndex], probFrequencies[zeroProbIndex]/(double)gridSize));
    }

    private void collectElevatedSampleBuckets()
    {
        mAllSampleBucketGroups = new HashMap();

        int totalCount = 0;
        double[][] probData = mBucketProbs.getData();

        for (int i = 0; i < mSampleCount; ++i)
        {
            List<Integer> bucketList = Lists.newArrayList();

            for (int j = 0; j < mBucketCount; ++j)
            {
                if (probData[j][i] > MAX_ELEVATED_PROB)
                {
                    continue;
                }

                bucketList.add(j);
                ++totalCount;
            }

            SampleData sample = mSampleData.get(i);
            sample.setElevatedBucketCounts(mElevatedCounts.getCol(i));

            if (!bucketList.isEmpty())
            {
                sample.setElevatedBuckets(bucketList);
                mAllSampleBucketGroups.put(i, bucketList);
            }
        }

        LOGGER.debug(String.format("samples with elevated buckets: count(%d perc=%.2f), buckets(%d perc=%.3f)",
                mAllSampleBucketGroups.size(), mAllSampleBucketGroups.size()/(double)mSampleCount,
                totalCount, totalCount/(double)(mBucketCount*mSampleCount)));
    }

    private void allocateSampleToGroup(int sampleId, BucketGroup bucketGroup, final double[] sampleCounts, boolean takeGroupRatios)
    {
        SampleData sample = mSampleData.get(sampleId);
        sample.addElevBucketGroup(bucketGroup);

        double allocCount = 0;
        double countsTotal = sumVector(sampleCounts);

        if(!takeGroupRatios)
        {
            allocCount = countsTotal;

            bucketGroup.addSample(sampleId, sampleCounts);
            sample.allocateBucketCounts(sampleCounts);
        }
        else
        {
            double[] allocCounts = getBestCountsAllocation(sampleId, bucketGroup.getBucketRatios(), sampleCounts);
            allocCount = sumVector(allocCounts);

            bucketGroup.addSample(sampleId, allocCounts);
            sample.allocateBucketCounts(allocCounts);
        }

//        LOGGER.debug(String.format("sample(%d) added to bg(%d) count(%s of %s perc=%.2f totalEvalPerc=%.2f)",
//                sampleId, bucketGroup.getId(), sizeToStr(allocCount), sizeToStr(countsTotal), allocCount/sample.getTotalCount(), sample.getAllocPercent()));
    }

    private double[] getBestCountsAllocation(int sampleId, final double[] bucketRatios, final double[] sampleCounts)
    {
        // must cap at the actual sample counts, but can go as high as the elevated probability allows
        final double[][] prData = mPermittedElevRange.getData();

        double minAlloc = 0;

        for(int i = 0; i < mBucketCount; ++i)
        {
            if(sampleCounts[i] == 0 || bucketRatios[i] == 0 || prData[i][sampleId] == 0)
                continue;

            double rangeHigh = sampleCounts[i] + prData[i][sampleId];
            double alloc = rangeHigh / bucketRatios[i];

            if(minAlloc == 0 || alloc < minAlloc)
                minAlloc = alloc;
        }

        double[] allocCounts = new double[mBucketCount];
        for(int i = 0; i < mBucketCount; ++i)
        {
            allocCounts[i] = minAlloc * bucketRatios[i];
        }

        return allocCounts;
    }

    private static int MIN_BUCKET_OVERLAP = 2;

    private void formExactBucketGroups()
    {
        // forms groups out of samples with similarly elevated buckets
        // returns true if groups were created or added to
        int groupsAdjusted = 0;
        int groupsCreated = 0;

        for (int samIndex1 = 0; samIndex1 < mSampleCount; ++samIndex1)
        {
            SampleData sample1 = mSampleData.get(samIndex1);

            final List<Integer> bl1 = sample1.getElevatedBuckets();

            if (bl1.isEmpty())
                continue;

            boolean addedToGroup = false;

            for (BucketGroup bucketGroup : mBucketGroups)
            {
                if(bucketGroup.hasSample(samIndex1))
                    continue;

                final List<Integer> commonBuckets = getMatchingList(bucketGroup.getBucketIds(), bl1);
                int commonBucketCount = commonBuckets.size();

                double groupMatch = commonBucketCount / (double)bucketGroup.getBucketIds().size();
                double sampleMatch = commonBucketCount / (double)bl1.size();
                double minMatch = min(groupMatch, sampleMatch);

                if (minMatch >= mHighRequiredMatch)
                {
                    List<Integer> combinedBuckets = getCombinedList(bl1, bucketGroup.getBucketIds());

                    double[] sc1 = extractBucketCountSubset(samIndex1, combinedBuckets);
                    double bcCss = calcSharedCSS(sc1, bucketGroup.getBucketCounts());

                    boolean withinProb = isWithinProbability(samIndex1, bucketGroup.getBucketRatios(), sc1, true);
                    boolean cssOK = bcCss >= mHighCssThreshold;

                    if(!cssOK || !withinProb)
                        continue;

                    /*
                    if(cssOK && !withinProb)
                    {
                        int[] result = checkPermittedSampleCounts(samIndex1, bucketGroup.getBucketRatios(), sc1, true, false);

                        double lowPerc = result[PI_LOW] / (double)result[PI_COUNT];
                        if(lowPerc >= 0.05)
                        {
                            double sampleBucketCount = sumVector(sc1);
                            LOGGER.debug(String.format("bg(%d) sample(%d) buckets(%d) elevCount(%.0f) has css(%.4f) prob(low=%d high=%d) mismatch",
                                    bucketGroup.getId(), samIndex1, commonBucketCount, sampleBucketCount, bcCss, result[PI_LOW], result[PI_HIGH]));
                            continue;
                        }
                    }
                    */

                    // group's buckets remain the same, neither increased nor refined, just add this new sample's counts
                    allocateSampleToGroup(samIndex1, bucketGroup, sc1, false);
                    ++groupsAdjusted;

                    LOGGER.debug(String.format("bg(%d) added sample(%d) with buckets(grp=%d sam=%d match=%d) css(%.4f prob=%s) totalSamples(%d)",
                            bucketGroup.getId(), samIndex1, bucketGroup.getBucketIds().size(), bl1.size(),
                            commonBuckets.size(), bcCss, withinProb ? "pass" : "fail", bucketGroup.getSampleIds().size()));

                    addedToGroup = true;
                        // break;
                }
            }

            if (addedToGroup)
            {
                // if this sample is now in a group, don't search amongst the other samples since they will each be test against the groups in turn
                continue;
            }

            for (int samIndex2 = samIndex1 + 1; samIndex2 < mSampleCount; ++samIndex2)
            {
                SampleData sample2 = mSampleData.get(samIndex2);

                final List<Integer> bl2 = sample2.getElevatedBuckets();

                if (bl2.isEmpty())
                    continue;

                if(mSampleWatchList.contains(samIndex1) && mSampleWatchList.contains(samIndex2))
                {
                    LOGGER.debug("specific sample");
                }

                final List<Integer> commonBuckets = getMatchingList(bl1, bl2);

                if(commonBuckets.size() < 2)
                    continue;

                int commonBucketCount = commonBuckets.size();
                double sample1Match = commonBucketCount / (double)bl1.size();
                double sample2Match = commonBucketCount / (double)bl2.size();
                double minMatch = min(sample1Match, sample2Match);

                if (minMatch >= mHighRequiredMatch && commonBucketCount >= MIN_BUCKET_OVERLAP)
                {
                    List<Integer> combinedBuckets = getCombinedList(bl1, bl2);
                    double[] sc1 = extractBucketCountSubset(samIndex1, combinedBuckets);
                    double[] sc2 = extractBucketCountSubset(samIndex2, combinedBuckets);
                    double bcCss = calcSharedCSS(sc1, sc2);

                    if(bcCss >= mHighCssThreshold)
                    {
                        BucketGroup bucketGroup = new BucketGroup(mBucketGroups.size());
                        bucketGroup.addBuckets(commonBuckets);

                        allocateSampleToGroup(samIndex1, bucketGroup, sc1, false);
                        allocateSampleToGroup(samIndex2, bucketGroup, sc2, false);

                        LOGGER.debug(String.format("added bg(%d) samples(%d and %d) with buckets(s1=%d s2=%d match=%d) css(%.4f)",
                                bucketGroup.getId(), samIndex1, samIndex2, bl1.size(), bl2.size(), commonBucketCount, bcCss));

                        mBucketGroups.add(bucketGroup);
                        ++groupsCreated;

                        // stop searching against other samples since they'll be tested against the groups in turn
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
                groupsCreated, groupsAdjusted, mBucketGroups.size(), mSampleCount - (groupsCreated * 2) - groupsAdjusted);
    }

    private void formBucketGroupsFromSamplePairSubsets()
    {
        // create unique groups from any subset of 2 samples' buckets if they are a close match
        // samples aren't necessarily allocated, only groups are made for the time-being
        LOGGER.debug("forming bucket groups from sample-pair reduced buckets");

        final double[][] elevData = mElevatedCounts.getData();

        // double startCssThreshold = 1 - (1 - mHighCssThreshold) * 0.5;
        double minCountThreshold = 20000; // mTotalCount * 0.0005; // was 0.00025, or about 10K
        int groupsAdded = 0;
        double[] combinedCounts = new double[mBucketCount];

        List<BucketGroup> bgList = Lists.newArrayList();

        for (int samIndex1 = 0; samIndex1 < mSampleCount; ++samIndex1)
        {
            SampleData sample1 = mSampleData.get(samIndex1);

            final List<Integer> bl1 = sample1.getElevatedBuckets();

            if (bl1.isEmpty())
                continue;

            for (int samIndex2 = samIndex1 + 1; samIndex2 < mSampleCount; ++samIndex2)
            {
                SampleData sample2 = mSampleData.get(samIndex2);

                final List<Integer> bl2 = sample2.getElevatedBuckets();

                if (bl2.isEmpty())
                    continue;

                if(mSampleWatchList.contains(samIndex1) && mSampleWatchList.contains(samIndex2))
                {
                    LOGGER.debug("specific sample");
                }

                List<Integer> commonBuckets = getMatchingList(bl1, bl2);
                int commonBucketCount = commonBuckets.size();

                if (commonBucketCount < 2)
                    continue;

                double[] sc1 = extractBucketCountSubset(samIndex1, commonBuckets);
                double[] sc2 = extractBucketCountSubset(samIndex2, commonBuckets);

                /*
                // check whether the overlapping buckets are covered by any of the existing groups
                // don't care about adding the samples, just want to know if the group exists
                BucketGroup sam1MatchesGroup = findMatchingBucketGroup(bgList, samIndex1, sc1, commonBuckets);
                BucketGroup sam2MatchesGroup = findMatchingBucketGroup(bgList, samIndex2, sc2, commonBuckets);

                if(sam1MatchesGroup != null && sam2MatchesGroup != null)
                    continue;
                */

                double allBucketsTotal = sample1.getElevatedCount() + sample2.getElevatedCount();
                double elevatedTotal = sumVector(sc1) + sumVector(sc2);

                // only create groups form overlapping buckets which contribute at least half the samples' mutation load
                if(elevatedTotal < allBucketsTotal * 0.5 || allBucketsTotal < minCountThreshold)
                    continue;

                double bcCss = calcSharedCSS(sc1, sc2);

                boolean addGroup = false;

                List<Integer> removedBuckets = Lists.newArrayList();

                if (bcCss >= mHighCssThreshold)
                {
                    addGroup = true;
                }
                else if (commonBucketCount > MIN_BUCKET_COUNT_OVERLAP)
                {
                    // attempt to find a match using less overlapping buckets
                    double cssThreshold = mHighCssThreshold;
                    // double cssThreshold = 1 - (1-startCssThreshold) * pow(bucketOverlapPerc, 2);

                    double[] cssResults = new double[commonBuckets.size()];

                    for (int i = 0; i < commonBuckets.size(); ++i)
                    {
                        Integer testBucket = commonBuckets.get(i);

                        sc1[testBucket] = 0;
                        sc2[testBucket] = 0;

                        // run CSS on this reduced set of buckets
                        cssResults[i] = calcSharedCSS(sc1, sc2);

                        // add the bucket back in ahead of the next test
                        sc1[testBucket] = elevData[testBucket][samIndex1];
                        sc2[testBucket] = elevData[testBucket][samIndex2];
                    }

                    List<Integer> sortedCssIndices = getSortedVectorIndices(cssResults, false);

                    // now remove each buckets one by one, with the ones having the largest effect to raise CSS first
                    for (int i = 0; i < sortedCssIndices.size(); ++i)
                    {
                        int commonBucketIndex = sortedCssIndices.get(i);
                        int testBucket = commonBuckets.get(commonBucketIndex);

                        sc1[testBucket] = 0;
                        sc2[testBucket] = 0;

                        elevatedTotal = sumVector(sc1) + sumVector(sc2);

                        // only consider groups from overlapping buckets which contribute at least half the samples' mutation load
                        if(elevatedTotal < allBucketsTotal * 0.5)
                            break;

                        // run CSS on this reduced set of buckets
                        bcCss = calcSharedCSS(sc1, sc2);

                        removedBuckets.add(testBucket);

                        if (bcCss >= cssThreshold)
                        {
                            addGroup = true;
                            break;
                        }

                        if(commonBucketCount - removedBuckets.size() <= MIN_BUCKET_COUNT_OVERLAP)
                            break;
                    }
                }

                if(!addGroup)
                    continue;

                List<Integer> matchedBuckets = Lists.newArrayList();

                if(removedBuckets.isEmpty())
                {
                    matchedBuckets.addAll(commonBuckets);
                }
                else
                {
                    for (Integer bucket : commonBuckets)
                    {
                        if (removedBuckets.contains(bucket))
                            continue;

                        matchedBuckets.add(bucket);
                    }
                }

                // with the reduced set of buckets, check again that a group doesn't already exist matching this new proposed one
                for(int i = 0; i < combinedCounts.length; ++i)
                {
                    combinedCounts[i] = sc1[i] + sc2[i];
                }

                if(findMatchingBucketGroup(bgList, samIndex1, combinedCounts, matchedBuckets) != null)
                {
                    addGroup = false;
                    continue;
                }

                BucketGroup bucketGroup = new BucketGroup(mBucketGroups.size() + bgList.size());

                bucketGroup.addBuckets(matchedBuckets);

                bucketGroup.addSample(samIndex1, sc1);
                bucketGroup.addSample(samIndex2, sc2);

                LOGGER.debug(String.format("added subset bg(%d) samples(%d and %d) with buckets(s1=%d and s2=%d common=%d matched=%d) counts(match=%s total=%s) css(%.4f)",
                        bucketGroup.getId(), samIndex1, samIndex2, bl1.size(), bl2.size(), commonBuckets.size(), matchedBuckets.size(),
                        sizeToStr(elevatedTotal), sizeToStr(allBucketsTotal), bcCss));

                bgList.add(bucketGroup);
                ++groupsAdded;
            }
        }

        mBucketGroups.addAll(bgList);

        LOGGER.debug("bucket groups created({}) total({})", groupsAdded, mBucketGroups.size());
    }

    private BucketGroup findMatchingBucketGroup(List<BucketGroup> bgList, int sampleId, double[] sampleCounts, List<Integer> bucketsToMatch)
    {
        if(bucketsToMatch.isEmpty())
            return null;

        // allow the group to have 25% more buckets but the CSS test is over the super set, so the additional ones can't be important
        double reqBucketOverlapPerc = 0.8;

        // match if the group has all of the required set and isn't much larger itself
        for (BucketGroup bucketGroup : bgList)
        {
            List<Integer> groupMatched = getMatchingList(bucketsToMatch, bucketGroup.getBucketIds());

            // all required buckets must be in the group
            if(groupMatched.size() != bucketsToMatch.size())
                continue;

            // the group must be primarily covered by this required set
            double bucketPerc = groupMatched.size() / (double)bucketGroup.getBucketIds().size();

            if (bucketPerc < reqBucketOverlapPerc)
                continue;

            // double cssThreshold = 1 - (1-startCssThreshold) * pow(bucketOverlapPerc, 2);

            // check CSS and if all good, then allocate these samples
            double groupCss = calcCSS(sampleCounts, bucketGroup.getBucketCounts(), false);
            if (groupCss >= mHighCssThreshold)
            {
                return bucketGroup;
            }
        }

        return null;
    }

    private void allocateSamplesToBucketGroups()
    {
        LOGGER.debug("allocating samples to {} bucket groups", mBucketGroups.size());

        // first score each group against the sample, then allocate mutuually exclusively to the highest scoring groups

        // scoring rules:
        // cumulative elevated counts across matching buckets

        // allocation rules:
        // skip a sample if it already fully allocated (99% or higher) or has no elevated buckets
        // allocate to a group if its bucket's counts equal or exceed the group's counts (for matching buckets)

        double fullyAllocatedPerc = 0.99;
        double minGroupAllocationPerc = 0.01;
        int maxListSize = 100;

        int GROUP_ID = 0;
        int ALLOC_ID = 1;

        double[] sampleCountAllocations = new double[mBucketCount];
        int SAMPLE_ID = 0;
        int UNALLOC_AMT = 1;

        int maxWorstSamples = 200;
        List<int[]> worstAllocatedSamples = Lists.newArrayList();

        // first clear all existing allocations
        for(BucketGroup bucketGroup : mBucketGroups)
        {
            bucketGroup.clearSamples();
        }

        for (int sampleId = 0; sampleId < mSampleCount; ++sampleId)
        {
            SampleData sample = mSampleData.get(sampleId);

            final List<Integer> samBuckets = sample.getElevatedBuckets();

            if(samBuckets.isEmpty())
                continue;

            if(sample.getAllocPercent() >= fullyAllocatedPerc)
                continue;

            if(mSampleWatchList.contains(sampleId))
            {
                // LOGGER.debug("specific sample");
            }

            List<int[]> bgScoreAllocations = Lists.newArrayList(); // a pair of BG index and allocation score

            double allocationThreshold = sample.getElevatedCount() * minGroupAllocationPerc;

            for (int grpIndex = 0; grpIndex < mBucketGroups.size(); ++grpIndex)
            {
                BucketGroup bucketGroup = mBucketGroups.get(grpIndex);

                final List<Integer> groupBuckets = bucketGroup.getBucketIds();
                final List<Integer> commonBuckets = getMatchingList(groupBuckets, samBuckets);

                if(mSampleWatchList.contains(bucketGroup.getId()))
                {
                    // LOGGER.debug("specific BG");
                }

                // ensure sample's buckets are covered by the group, even if the sample has more
                double groupPerc = commonBuckets.size() / (double) groupBuckets.size();
                if(groupPerc < mHighRequiredMatch)
                    continue;

                double[] sampleCounts = extractBucketCountSubset(sampleId, groupBuckets);

                double[] allocCounts = getBestCountsAllocation(sampleId, bucketGroup.getBucketRatios(), sampleCounts);
                double allocCountTotal = sumVector(allocCounts);

                if(allocCountTotal < allocationThreshold)
                    continue;

                if(allocCountTotal > sample.getElevatedCount() * 1.2)
                {
                    LOGGER.warn(String.format("sample(%d) allocation(%.0f) exceeds elevatedTotal(%.0f)", sampleId, allocCountTotal, sample.getElevatedCount()));
                }

                int[] bgAndScore = { grpIndex, (int)round(allocCountTotal) };

                int index = 0;
                for(; index < bgScoreAllocations.size(); ++index)
                {
                    if(allocCountTotal > bgScoreAllocations.get(index)[ALLOC_ID])
                        break;

                    if(index >= maxListSize)
                        break;
                }

                if(index < maxListSize)
                {
                    bgScoreAllocations.add(index, bgAndScore);
                }
            }

            if(bgScoreAllocations.isEmpty())
            {
                LOGGER.debug(String.format("sample(%d) found no matching bucket groups, elev buckets(%d) count(%s)",
                        sampleId, sample.getElevatedBuckets().size(), sizeToStr(sample.getElevatedCount())));
                continue;
            }

            LOGGER.debug("sample({}) has {} non-zero possible bucket groups", sampleId, bgScoreAllocations.size());

            int iterations = 0;
            int maxGroups = 10;

            for(int[] bgAndScore : bgScoreAllocations)
            {
                // get the next highest group allocation
                int grpIndex = bgAndScore[GROUP_ID];
                double allocCounts = bgAndScore[ALLOC_ID];

                if(allocCounts == 0)
                    break;

                BucketGroup bucketGroup = mBucketGroups.get(grpIndex);

                final List<Integer> groupBuckets = bucketGroup.getBucketIds();
                double[] bucketRatios = bucketGroup.getBucketRatios();

                // apply to the remaining unallocated elevated counts for this sample
                for(int j = 0; j < mBucketCount; ++j)
                {
                    sampleCountAllocations[j] = 0;

                    if(!groupBuckets.contains(j))
                        continue;

                    sampleCountAllocations[j] = allocCounts * bucketRatios[j];
                }

                double prevAllocPerc = sample.getAllocPercent();
                sample.allocateBucketCounts(sampleCountAllocations);

                if(sample.getAllocPercent() >= prevAllocPerc + minGroupAllocationPerc)
                {
                    // note: no attempt to remove tiny allocations from the sample, since not interested in actual residuals
                    bucketGroup.addSample(sampleId, sampleCountAllocations);
                    sample.addElevBucketGroup(bucketGroup);

                    LOGGER.debug(String.format("sample(%d) added to bg(%d buckets=%d sam=%d) count(%s of %s) allocatedPerc(%.3f -> %.3f)",
                            sampleId, bucketGroup.getId(), groupBuckets.size(), sample.getElevatedBuckets().size(),
                            sizeToStr(allocCounts), sizeToStr(sample.getElevatedCount()), prevAllocPerc, sample.getAllocPercent()));

                    if (sample.getAllocPercent() >= fullyAllocatedPerc)
                        break;

                    if(sample.getElevBucketGroups().size() >= maxGroups)
                        break;
                }

                ++iterations;
            }

            LOGGER.debug(String.format("sample(%d) allocated to groups(%d), elevated buckets(%d) count(%.3f of %s) iter(%d)",
                    sampleId, sample.getElevBucketGroups().size(), sample.getElevatedBuckets().size(),
                    sample.getAllocPercent(), sizeToStr(sample.getElevatedCount()), iterations));

            int unallocTotal = (int)round((1 - sample.getAllocPercent()) * sample.getElevatedCount());
            int worstIndex = 0;
            for(; worstIndex < worstAllocatedSamples.size(); ++worstIndex)
            {
                if(unallocTotal > worstAllocatedSamples.get(worstIndex)[UNALLOC_AMT])
                    break;
            }

            if(worstIndex <= maxWorstSamples)
            {
                int[] worstData = { sampleId, unallocTotal };
                worstAllocatedSamples.add(worstIndex, worstData);
            }
        }

        // log worst x samples
        for(int worstIndex = 0; worstIndex < min(worstAllocatedSamples.size(), maxWorstSamples); ++worstIndex)
        {
            final int[] worstData = worstAllocatedSamples.get(worstIndex);
            final SampleData sample = mSampleData.get(worstData[SAMPLE_ID]);

            LOGGER.debug(String.format("%d: worst sample(%d: %s) cancer(%s) unallocated(%s of %s, perc=%.3f)",
                    worstIndex, sample.Id, sample.getSampleName(), sample.getCancerType(),
                    sizeToStr(worstData[UNALLOC_AMT]), sizeToStr(sample.getElevatedCount()), (1 - sample.getAllocPercent())));
        }
    }

    private void formBucketGroupsFromSubsets()
    {
        LOGGER.debug("forming bucket groups from like-pairs of samples with reduced overlap");

        final double[][] elevData = mElevatedCounts.getData();

        double startCssThreshold = 1 - (1 - mHighCssThreshold) * 0.5;
        double minCountThreshold = mTotalCount * 0.00025;
        int groupsAdded = 0;

        for (int samIndex1 = 0; samIndex1 < mSampleCount; ++samIndex1)
        {
            final List<Integer> bl1 = mAllSampleBucketGroups.get(samIndex1);

            if (bl1 == null)
                continue;

            List<Integer> bgIndices = getSampleBucketGroupIndices(samIndex1);

            List<BucketGroup> sam1Groups = Lists.newArrayList();
            for (Integer bgIndex : bgIndices)
            {
                sam1Groups.add(mBucketGroups.get(bgIndex));
            }

            double sam1Total = sumVector(mElevatedCounts.getCol(samIndex1));

            for (int samIndex2 = samIndex1 + 1; samIndex2 < mSampleCount; ++samIndex2)
            {
                final List<Integer> bl2 = mAllSampleBucketGroups.get(samIndex2);

                if (bl2 == null)
                    continue;

                if(mSampleWatchList.contains(samIndex1) && mSampleWatchList.contains(samIndex2))
                {
                    LOGGER.debug("specific sample");
                }

                List<Integer> commonBuckets = getMatchingList(bl1, bl2);
                int commonBucketCount = commonBuckets.size();

                if (commonBucketCount < 2)
                    continue;

                // if (commonBucketCount < MIN_BUCKET_COUNT_OVERLAP)
                //      continue;

                double[] sc1 = extractBucketCountSubset(samIndex1, commonBuckets);
                double[] sc2 = extractBucketCountSubset(samIndex2, commonBuckets);

                // check whether the overlapping buckets are covered by any of the first sample's groups
                boolean sam1GroupsContainsList = false;
                boolean sam2AddedToGroup = false;
                for(final BucketGroup bucketGroup : sam1Groups)
                {
                    if(bucketGroup.hasBuckets(commonBuckets))
                    {
                        double groupCss = calcSharedCSS(sc2, bucketGroup.getBucketCounts());
                        if(groupCss >= mHighCssThreshold && isWithinProbability(samIndex2, bucketGroup.getBucketRatios(), sc2, true))
                        {
                            bucketGroup.addSample(samIndex2, sc2);

                            LOGGER.debug(String.format("bg(%d) added sample(%d) with buckets(grp=%d sam=%d match=%d) css(%.4f) totalSamples(%d)",
                                    bucketGroup.getId(), samIndex2, bucketGroup.getBucketIds().size(), bl2.size(), commonBuckets.size(),
                                    groupCss, bucketGroup.getSampleIds().size()));

                            sam2AddedToGroup = true;
                        }

                        sam1GroupsContainsList = true;
                        break;
                    }
                }

                if(sam2AddedToGroup)
                    continue;

                List<Integer> bg2Indices = getSampleBucketGroupIndices(samIndex2);

                boolean sam2GroupsContainsList = false;
                for(Integer bgIndex : bg2Indices)
                {
                    final BucketGroup group = mBucketGroups.get(bgIndex);

                    if(group.hasBuckets(commonBuckets))
                    {
                        sam2GroupsContainsList = true;
                        break;
                    }
                }

//                if(sam2GroupsContainsList)
//                {
//                    LOGGER.debug("sample1({}) covered by one of sample2({})'s groups", samIndex1, samIndex2);
//                }

                if(sam1GroupsContainsList && sam2GroupsContainsList)
                     continue;

                double allBucketsTotal = sam1Total + sumVector(mElevatedCounts.getCol(samIndex2));
                double elevatedTotal = sumVector(sc1) + sumVector(sc2);

                // only create groups form overlapping buckets which contribute at least half the samples' mutation load
                if(elevatedTotal < allBucketsTotal * 0.5 || allBucketsTotal < minCountThreshold)
                    continue;

                double bcCss = calcSharedCSS(sc1, sc2);

                boolean addGroup = false;

                List<Integer> removedBuckets = Lists.newArrayList();

                if (bcCss >= startCssThreshold)
                {
                    addGroup = true;
                }
                else if (commonBucketCount > MIN_BUCKET_COUNT_OVERLAP)
                {
                    // attempt to find a match using less overlapping buckets
                    int currentBucketCount = commonBucketCount;

                    List<Integer> testBuckets = Lists.newArrayList();
                    testBuckets.addAll(commonBuckets);

                    while (!addGroup && currentBucketCount > MIN_BUCKET_COUNT_OVERLAP)
                    {
                        double bestCss = 0;
                        int bestTestBucket = 0;

                        --currentBucketCount;

                        double cssThreshold = mHighCssThreshold;
                        // double cssThreshold = 1 - (1-startCssThreshold) * pow(bucketOverlapPerc, 2);

                        for (Integer removedBucket : removedBuckets)
                        {
                            sc1[removedBucket] = 0;
                            sc2[removedBucket] = 0;
                        }

                        for (Integer testBucket : testBuckets)
                        {
//                            if(removedBuckets.contains(testBucket))
//                                continue;

                            sc1[testBucket] = 0;
                            sc2[testBucket] = 0;

                            // run CSS on this reduced set of buckets
                            bcCss = calcSharedCSS(sc1, sc2);

                            // add the bucket back in ahead of the next test
                            sc1[testBucket] = elevData[testBucket][samIndex1];
                            sc2[testBucket] = elevData[testBucket][samIndex2];

                            if (bcCss > bestCss)
                            {
                                bestCss = bcCss;
                                bestTestBucket = testBucket;
                            }

                            if (bcCss >= cssThreshold)
                                break;
                        }

                        removedBuckets.add(bestTestBucket);
                        sc1[bestTestBucket] = 0;
                        sc2[bestTestBucket] = 0;

                        for(int i = 0; i < testBuckets.size(); ++i)
                        {
                            if(testBuckets.get(i) == bestTestBucket)
                            {
                                testBuckets.remove(i);
                                break;
                            }
                        }

                        elevatedTotal = sumVector(sc1) + sumVector(sc2);

                        // only consider groups from overlapping buckets which contribute at least half the samples' mutation load
                        if(elevatedTotal < allBucketsTotal * 0.5)
                            break;

                        if (bestCss >= cssThreshold)
                        {
                            addGroup = true;
                            break;
                        }
                    }
                }

                if(addGroup)
                {
                    boolean addedToGroup = false;

                    List<Integer> matchedBuckets = Lists.newArrayList();
                    for (Integer bucket : commonBuckets)
                    {
                        if (removedBuckets.contains(bucket))
                            continue;

                        matchedBuckets.add(bucket);
                    }

                    // check that a group matching this hasn't just been added before creating a new one
                    double groupCss = 0;

                    for (BucketGroup bucketGroup : mBucketGroups)
                    {
                        if(bucketGroup.getBucketIds().size() != matchedBuckets.size() || !bucketGroup.hasBuckets(matchedBuckets))
                            continue;

                        // check CSS and if all good, then add these samples
                        if (!bucketGroup.hasSample(samIndex1))
                        {
                            groupCss = calcSharedCSS(sc1, bucketGroup.getBucketCounts());
                            if (groupCss >= mHighCssThreshold)
                            {
                                bucketGroup.addSample(samIndex1, sc1);

                                LOGGER.debug(String.format("bg(%d) added sample(%d) with buckets(grp=%d sam=%d match=%d) css(%.4f) totalSamples(%d)",
                                        bucketGroup.getId(), samIndex1, bucketGroup.getBucketIds().size(), bl1.size(), matchedBuckets.size(), groupCss, bucketGroup.getSampleIds().size()));

                                addedToGroup = true;
                            }
                        }

                        // check CSS and if all good, then add these samples
                        if (!bucketGroup.hasSample(samIndex2))
                        {
                            groupCss = calcSharedCSS(sc2, bucketGroup.getBucketCounts());
                            if (groupCss >= mHighCssThreshold)
                            {
                                bucketGroup.addSample(samIndex2, sc2);

                                LOGGER.debug(String.format("bg(%d) added sample(%d) with buckets(grp=%d sam=%d match=%d) css(%.4f) totalSamples(%d)",
                                        bucketGroup.getId(), samIndex2, bucketGroup.getBucketIds().size(), bl2.size(), matchedBuckets.size(), groupCss, bucketGroup.getSampleIds().size()));

                                addedToGroup = true;
                            }
                        }

                        if(addedToGroup)
                            break;
                    }

                    if(!addedToGroup)
                    {
                        BucketGroup bucketGroup = new BucketGroup(mBucketGroups.size());

                        bucketGroup.addBuckets(matchedBuckets);
                        bucketGroup.addSample(samIndex1, sc1);
                        bucketGroup.addSample(samIndex2, sc2);

                        LOGGER.debug(String.format("added subset bg(%d) samples(%d and %d) with buckets(s1=%d and s2=%d common=%d matched=%d) counts(match=%s total=%s) css(%.4f)",
                                bucketGroup.getId(), samIndex1, samIndex2, bl1.size(), bl2.size(), commonBuckets.size(), matchedBuckets.size(),
                                sizeToStr(elevatedTotal), sizeToStr(allBucketsTotal), bcCss));

                        mBucketGroups.add(bucketGroup);
                        sam1Groups.add(bucketGroup);
                        ++groupsAdded;
                    }
                }
            }
        }

        LOGGER.debug("bucket groups created({}) total({})", groupsAdded, mBucketGroups.size());
    }

    private void broadenBucketGroupDefinitions()
    {
        LOGGER.debug("broadening bucket groups with extra buckets");

        for (BucketGroup bucketGroup : mBucketGroups)
        {
            // starting with the current set of buckets, attempt to add in a new bucket as long as the CSS remains above the threshold
            // but only apply if the samples in the group show evidence of being elevated in the new buckets being considered
            final List<Integer> startBuckets = bucketGroup.getBucketIds();
            final List<Integer> sampleIds = bucketGroup.getSampleIds();
            List<Integer> workingBuckets = Lists.newArrayList();
            workingBuckets.addAll(startBuckets);

            // populate the matrix
            NmfMatrix sampleCounts = new NmfMatrix(mBucketCount, sampleIds.size());
            double[][] scData = sampleCounts.getData();
            final double[][] refData = mElevatedCounts.getData();
            final double[][] bmrData = mBucketMedianRatios.getData();

            for(int i = 0; i < mBucketCount; ++i)
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
            for (int testBucket = 0; testBucket < mBucketCount; ++testBucket)
            {
                if (workingBuckets.contains(testBucket))
                    continue;

                // check if this new bucket tested shows evidence of being elevated in the samples
                boolean areElevated = true;
                for (Integer sampleId : sampleIds)
                {
                    if(bmrData[testBucket][sampleId] < MIN_ELEVATION_RATIO)
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
                    lastTestedBucket = -1; // to retain this new bucket for future tests

                    // add this bucket
                    bucketGroup.addBucket(testBucket, sampleCounts.getRow(testBucket), false);
                }
            }

            if(bucketsAdded > 0)
            {
                LOGGER.debug("bg({}) samples({}) broadened buckets(init={} new={} added={}) ",
                        bucketGroup.getId(), bgSampleCount, bucketGroup.getBucketIds().size() - bucketsAdded,
                        bucketGroup.getBucketIds().size(), bucketsAdded);
            }
        }
    }

    private void checkAllSamplesAgainstBucketGroups()
    {
        LOGGER.debug("checking all samples against {} bucket groups", mBucketGroups.size());

        int samplesAdded = 0;

        final double[][] bmrData = mBucketMedianRatios.getData();

        for (BucketGroup bucketGroup : mBucketGroups)
        {
            if(bucketGroup.getBucketIds().size() == 1)
                continue; // can't calc CSS using a single bucket

            final List<Integer> groupBuckets = bucketGroup.getBucketIds();

            for (int sampleId = 0; sampleId < mSampleCount; ++sampleId)
            {
                if(bucketGroup.hasSample(sampleId))
                    continue; // ignore those already assigned

                if(mSampleWatchList.contains(sampleId))
                {
                    LOGGER.debug("specific sample");
                }

                double[] samCounts = extractBucketCountSubset(sampleId, groupBuckets);
                double bcCss = calcSharedCSS(samCounts, bucketGroup.getBucketCounts());

                boolean withinProb = isWithinProbability(sampleId, bucketGroup.getBucketRatios(), samCounts, true);

                boolean cssOK = bcCss >= mHighCssThreshold;

                if(!cssOK && !withinProb)
                    continue;

                final List<Integer> samBuckets = mAllSampleBucketGroups.get(sampleId);

                if(samBuckets == null)
                    continue; // must have at least 1 elevated bucket

                final List<Integer> commonBuckets = getMatchingList(groupBuckets, samBuckets);

                if(cssOK && !withinProb)
                {
                    int[] result = checkPermittedSampleCounts(sampleId, bucketGroup.getBucketRatios(), samCounts, true, false);

                    double lowPerc = result[PI_LOW] / (double)result[PI_COUNT];
                    if(lowPerc >= 0.05)
                    {
                        // double sampleBucketCount = sumVector(samCounts);
                        // LOGGER.debug(String.format("bg(%d) sample(%d) buckets(%d) elevCount(%.0f) has css(%.4f) prob(low=%d high=%d) mismatch",
                        //         bucketGroup.getId(), sampleId, bucketGroup.getBucketIds().size(), sampleBucketCount, bcCss, result[PI_LOW], result[PI_HIGH]));
                        continue;
                    }
                }

                double groupBucketOverlap = commonBuckets.size() / (double)groupBuckets.size();

                if(groupBucketOverlap < mHighRequiredMatch)
                    continue;

                final List<Integer> bucketsNotInSample = getDiffList(groupBuckets, samBuckets);

                if(!bucketsNotInSample.isEmpty())
                {
                    boolean allBucketsElevatedBmr = true;
                    for (Integer bucket : bucketsNotInSample)
                    {
                        if (bmrData[bucket][sampleId] < MIN_ELEVATION_RATIO)
                        {
                            allBucketsElevatedBmr = false;
                            break;
                        }
                    }

                    if (!allBucketsElevatedBmr)
                        continue;
                }

                // don't need to check buckets that are unique to the sample since it can belong to more than one group
                bucketGroup.addSample(sampleId, samCounts);
                ++samplesAdded;

                LOGGER.debug(String.format("bg(%d) adding sample(%d) with buckets(grp=%d samElev=%d grpMatch=%.2f bmr=pass) css(%.4f prob=%s) totalSamples(%d)",
                        bucketGroup.getId(), sampleId, groupBuckets.size(), samBuckets.size(), groupBucketOverlap,
                        bcCss, withinProb ? "pass" : "fail", bucketGroup.getSampleIds().size()));
            }
        }

        LOGGER.debug("added {} samples to existing bucket groups", samplesAdded);
    }

    private void mergeBucketGroups()
    {
        removeEmptyBucketGroups();

        LOGGER.debug("checking merge for {} bucket groups", mBucketGroups.size());

        double reqAllMatchPercent = 0.85; // bucket overlap including those from widening
        double reqInitMatchPercent = 0.95; // bucket overlap from initial set of elevated buckets

        // what bucket similarity should be considered to merge 2 bucket groups?
        // mustn't lose unique buckets, eg if not explained by some other bucket group
        // and the collapsing-sample routine to come will also do some of the same work

        // rules: merge smaller sample count into larger
        // rules: merge only if the diff buckets in the smaller group can be explained by another BG exactly
        // OR don't merge if the differing buckets are from the initial sets

        int bgIndex1 = 0;
        while (bgIndex1 < mBucketGroups.size())
        {
            BucketGroup bg1 = mBucketGroups.get(bgIndex1);

            if(bg1.getSampleIds().size() < 2 || bg1.getBucketIds().size() < 2)
            {
                mBucketGroups.remove(bgIndex1);
                continue;
            }

            boolean removeBg1 = false;

            int bgIndex2 = bgIndex1+1;
            while (bgIndex2 < mBucketGroups.size())
            {
                BucketGroup bg2 = mBucketGroups.get(bgIndex2);

                // first check CSS on the common buckets only
                double bcCss = calcSharedCSS(bg1.getBucketCounts(), bg2.getBucketCounts());

                final List<Integer> commonAllBuckets = getMatchingList(bg1.getBucketIds(), bg2.getBucketIds());
                boolean hasMatchingBuckets = (bg1.getBucketIds().size() == bg2.getBucketIds().size()) && (bg1.getBucketIds().size() == commonAllBuckets.size());

                // if buckets are identical, have a lower CSS check, consistent with how proposed sigs will be compared
                if ((hasMatchingBuckets && bcCss < mSigCssThreshold) || (!hasMatchingBuckets && bcCss < mHighCssThreshold))
                {
                    ++bgIndex2;
                    continue;
                }

                boolean removeBg2 = false;

                final List<Integer> commonInitBuckets = getMatchingList(bg1.getInitialBucketIds(), bg2.getInitialBucketIds());

                int commonInitCount = commonInitBuckets.size();
                int commonAllCount = commonAllBuckets.size();
                int bc1 = bg1.getBucketIds().size();
                int bc2 = bg2.getBucketIds().size();

                double minMatchInit = min(commonInitCount/(double)bc1, commonInitCount/(double)bc2);
                double minMatchAll = min(commonAllCount/(double)bc1, commonAllCount/(double)bc2);

                if (minMatchAll < reqAllMatchPercent || minMatchInit < reqInitMatchPercent)
                {
                    // check CSS again but this time on the super set of buckets and not ignoring zeros
                    bcCss = calcCSS(bg1.getBucketCounts(), bg2.getBucketCounts(), false);

                    if(bcCss < mHighCssThreshold)
                    {
                        ++bgIndex2;
                        continue;
                    }

                    LOGGER.debug(String.format("bg(%d) vs bg(%d) differing buckets(bg1=%d bg2=%d match=%d) but high css(%.4f) on superset",
                            bg1.getId(), bg2.getId(), bg1.getBucketIds().size(), bg2.getBucketIds().size(), commonAllCount, bcCss));

//                    if(minMatchAll >= reqAllMatchPercent * 0.8 && bcCss >= 0.98)
//                    {
//                        LOGGER.debug(String.format("bg(%d) vs bg(%d) close-to-merge buckets(bg1=%d bg2=%d match=%d) css(%.4f)",
//                                bg1.getId(), bg2.getId(), bg1.getBucketIds().size(), bg2.getBucketIds().size(), commmonBucketCount, bcCss));
//                    }
//                    else if(minMatch >= 0.8 && bcCss < 0.5)
//                    {
//                        // check for overlapping different processes?
//                        LOGGER.debug(String.format("bg(%d) vs bg(%d) very diff despite buckets(bg1=%d bg2=%d match=%d) css(%.4f)",
//                                bg1.getId(), bg2.getId(), bg1.getBucketIds().size(), bg2.getBucketIds().size(), commmonBucketCount, bcCss));
//                    }

                }

                int bg1SC = bg1.getSampleIds().size();
                int bg2SC = bg2.getSampleIds().size();

                if(minMatchAll < 1 || minMatchInit < 1)
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
                        bg1.reduceToBucketSet(commonAllBuckets);
                        removeBg2 = true;
                    }
                }
                else
                {
                    // both groups have the same buckets so merging either direction is fine
                    bg1.merge(bg2.getSampleIds(), bg2.getBucketCounts());
                    removeBg2 = true;
                }

                LOGGER.debug(String.format("bg(%d) %s bg(%d) buckets(bg1=%d bg2=%d matchInit=%.2f matchAll=%.2f) css(%.4f) samples(bg1=%d bg2=%d total=%d)",
                        bg1.getId(), removeBg2 ? "merges in" : "merged into", bg2.getId(), bg1.getBucketIds().size(), bg2.getBucketIds().size(),
                        minMatchInit, minMatchAll, bcCss, bg1SC, bg2SC, bg1SC + bg2SC));

                if(removeBg2)
                {
                    mBucketGroups.remove(bgIndex2);
                    continue;
                }
                else
                {
                    // BG 1 will be removed
                    break;
                }
            }

            if(removeBg1)
                mBucketGroups.remove(bgIndex1);
            else
                ++bgIndex1;
        }

        LOGGER.debug("post-merge has {} bucket groups", mBucketGroups.size());
    }

    private void collapseBucketSamples()
    {
        LOGGER.debug("collapsing samples between overlapping bucket group");

        // if a sample is allocated to 2 bucket groups, where one is wholy covered by another,
        // then remove it from the smaller group
        // this will then lead to some groups being empty and then removed
        for (int bgIndex1 = 0; bgIndex1 < mBucketGroups.size(); ++bgIndex1)
        {
            BucketGroup bg1 = mBucketGroups.get(bgIndex1);

            for (int bgIndex2 = bgIndex1 + 1; bgIndex2 < mBucketGroups.size(); ++bgIndex2)
            {
                BucketGroup bg2 = mBucketGroups.get(bgIndex2);

                final List<Integer> commonBuckets = getMatchingList(bg1.getBucketIds(), bg2.getBucketIds());
                int commmonBucketCount = commonBuckets.size();

                if (commmonBucketCount != bg1.getBucketIds().size() && commmonBucketCount != bg2.getBucketIds().size())
                    continue;

                boolean bg1IsSubset = commmonBucketCount == bg1.getBucketIds().size();

                double bcCss = calcSharedCSS(bg1.getBucketCounts(), bg2.getBucketCounts());

                // if buckets are identical, have a lower CSS check, consistent with how proposed sigs will be compared
                if (bcCss < mHighCssThreshold)
                {
                    if(bcCss >= mSigCssThreshold)
                    {
//                        LOGGER.debug(String.format("bg(%d) %s with bg(%d) buckets(bg1=%d bg2=%d match=%d) css(%.4f) samples(bg1=%d bg2=%d) too low to collapse",
//                                bg1.getId(), bg1IsSubset ? "subset" : "superset", bg2.getId(), bg1.getBucketIds().size(), bg2.getBucketIds().size(),
//                                commmonBucketCount, bcCss, bg1.getSampleIds().size(), bg2.getSampleIds().size()));
                    }
                    continue;
                }

                // move samples out of this bucket group if they're in the larger group

                final List<Integer> mainSamples = bg1IsSubset ? bg2.getSampleIds() : bg1.getSampleIds();
                List<Integer> subsetSamples = bg1IsSubset ? bg1.getSampleIds() : bg2.getSampleIds();

                int samplesRemoved = 0;
                int samIndex = 0;
                while(samIndex < subsetSamples.size())
                {
                    Integer sampleId = subsetSamples.get(samIndex);

                    if(mainSamples.contains(sampleId))
                    {
                        subsetSamples.remove(samIndex);
                        ++samplesRemoved;
                    }
                    else
                    {
                        ++samIndex;
                    }
                }

                if(samplesRemoved > 0)
                {
                    LOGGER.debug(String.format("bg(%d) %s with bg(%d) buckets(bg1=%d bg2=%d match=%d) css(%.4f) samples(bg1=%d bg2=%d switched=%d)",
                            bg1.getId(), bg1IsSubset ? "subset" : "superset", bg2.getId(), bg1.getBucketIds().size(), bg2.getBucketIds().size(),
                            commmonBucketCount, bcCss, bg1.getSampleIds().size(), bg2.getSampleIds().size(), samplesRemoved));
                }
            }
        }

        removeEmptyBucketGroups();
    }

    private void removeEmptyBucketGroups()
    {
        // finally remove any bucket groups which now have too few samples or buckets
        int bgIndex = 0;
        int bgRemoved = 0;
        while(bgIndex < mBucketGroups.size())
        {
            BucketGroup bucketGroup = mBucketGroups.get(bgIndex);

            if(bucketGroup.getSampleIds().size() < 2 || bucketGroup.getBucketIds().size() < 2)
            {
                mBucketGroups.remove(bgIndex);
                ++bgRemoved;
            }
            else
            {
                ++bgIndex;
            }
        }

        if(bgRemoved > 0)
        {
            LOGGER.debug("removed {} bucket groups with < 2 samples or < 2 buckets", bgRemoved);
        }
    }

    private void recalcBucketData()
    {
        // post merging, recalc bucket counts across all samples for a group
        // and work out purity as percentage of all sample buckets that are elevated
        final double[][] elevData = mElevatedCounts.getData();
        final double[][] bmrData = mBucketMedianRatios.getData();

        int totalGridSize = mBucketCount * mSampleCount;
        double avgCount = mTotalCount/totalGridSize;

        for (final BucketGroup bucketGroup : mBucketGroups)
        {
            final List<Integer> sampleIds = bucketGroup.getSampleIds();
            final List<Integer> bucketIds = bucketGroup.getBucketIds();
            List<Double> sampleCountTotals = Lists.newArrayList();

            // recalc the bucket count totals since merging did not attempt to handle double-counting sample's counts
            double[] bucketCounts = new double[mBucketCount];

            int lowProbCount = 0;
            double totalBmr = 0;

            for (Integer sampleId : sampleIds)
            {
                SampleData sample = mSampleData.get(sampleId);

                final List<Integer> elevatedBuckets = sample.getElevatedBuckets();

                double sampleMutLoad = 0; // just for the applicable buckets

                for (Integer bucketId : bucketIds)
                {
                    double sbCount = elevData[bucketId][sampleId];
                    sampleMutLoad += sbCount;
                    bucketCounts[bucketId] += sbCount;

                    if (elevatedBuckets != null && elevatedBuckets.contains(bucketId))
                    {
                        ++lowProbCount;
                    }

                    totalBmr += bmrData[bucketId][sampleId];
                }

                sampleCountTotals.add(sampleMutLoad);
            }

            bucketGroup.setBucketCounts(bucketCounts);
            bucketGroup.setSampleCountTotals(sampleCountTotals);

            int groupSize = bucketGroup.getSize();
            double avgItemCount = sumVector(bucketCounts) / groupSize;
            double avgBmr = totalBmr / groupSize;
            double loadFactor = max(log(avgBmr * avgItemCount/avgCount), 0.1);

            double purity = lowProbCount / (double)groupSize;
            bucketGroup.setPurity(purity);
            bucketGroup.setLoadFactor(loadFactor);

            LOGGER.debug(String.format("bg(%d) size(%d samples=%d buckets=%d) purity(%.2f) score(%.0f) loadFactor(%.1f) avgItemCount(%.0f vs avg=%.0f) avgBmr(%.1f)",
                    bucketGroup.getId(), bucketGroup.getSize(), sampleIds.size(), bucketIds.size(),
                    purity, bucketGroup.calcScore(), loadFactor, avgItemCount, avgCount, avgBmr));
        }
    }

    private void analyseGroupsVsExtData(boolean useElevated, boolean verbose)
    {
        if(mExtSampleData == null)
            return;

        List<BucketGroup> bgList = useElevated ? mBucketGroups : mBackgroundGroups;

        LOGGER.debug("analysing {} bucket groups", bgList.size());

        // check counts for each external data category against the samples in each bucket group

        // want to know if:
        // a) the samples have 1 or more identical categoroes and
        // b) the proportion of these categories out of the cohort (ie missing other samples, linked to other bucket groups)

        for (BucketGroup bucketGroup : bgList)
        {
            final List<Integer> sampleIds = bucketGroup.getSampleIds();
            int sampleCount = sampleIds.size();

            // clear any previous
            bucketGroup.setCancerType("");
            bucketGroup.setEffects("");

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

                if(verbose)
                {
                    LOGGER.debug(String.format("bg(%d) category(%s) count(%d of %d) perc(group=%.3f category=%.3f) buckets(%d)",
                            bucketGroup.getId(), catName, catSampleCount, sampleCount, samplesPerc, catPerc, bucketGroup.getBucketIds().size()));
                }

                if(extractCategoryName(catName).equals(CATEGORY_CANCER_TYPE))
                {
                    String cancerType = extractCategoryValue(catName);

                    if(samplesPerc >= DOMINANT_CATEGORY_PERCENT)
                    {
                        bucketGroup.setCancerType(String.format("%s=%.2f", cancerType, samplesPerc));
                    }
                    else
                    {
                        String cancerTypeStr = bucketGroup.getCancerType();

                        if(!cancerTypeStr.isEmpty())
                            cancerTypeStr += ";";

                        cancerTypeStr += String.format("%s=%.2f", cancerType, samplesPerc);

                        bucketGroup.setCancerType(cancerTypeStr);
                    }
                }
                else
                {
                    String effects = bucketGroup.getEffects();

                    if(!effects.isEmpty())
                        effects += ";";

                    effects += String.format("%s=%.2f", catName, samplesPerc);

                    bucketGroup.setEffects(effects);
                }
            }
        }
    }

    private void logBucketGroups(boolean useElevated, boolean verbose)
    {
        List<BucketGroup> bgList = useElevated ? mBucketGroups : mBackgroundGroups;

        // sort by score before logging
        Collections.sort(bgList);

        // log top groups
        int maxToLog = verbose ? bgList.size() : min(bgList.size(), 40);
        int minScore = verbose ? 0 : 10; // 4 samples x 5 buckets, or 2 x 10, and adjusted for purity

        if(verbose)
        {
            LOGGER.debug("logging all bucket groups of total({})", bgList.size());
        }
        else
        {
            LOGGER.debug("logging top {} bucket groups of total({})", maxToLog, bgList.size());
        }

        double totalCount = useElevated ? mElevatedCount : mBackgroundCounts.sum();
        final String groupType = useElevated ? "elevated" : "background";

        for (int i = 0; i < maxToLog; ++i)
        {
            BucketGroup bucketGroup = bgList.get(i);

            if (bucketGroup.calcScore() < minScore)
                break;

            if (verbose)
            {
                String bucketIdsStr = "";

                if(useElevated)
                {
                    // log top 20 initial buckets, then any extra ones added through the broadending routine
                    final List<Integer> initBuckets = bucketGroup.getInitialBucketIds();
                    List<Integer> descBucketRatioIndices = getSortedVectorIndices(bucketGroup.getBucketRatios(), false);
                    for(int j = 0; j < min(20, descBucketRatioIndices.size()); ++j)
                    {
                        Integer bucket = descBucketRatioIndices.get(j);

                        if(!initBuckets.contains(bucket))
                            continue;

                        if (!bucketIdsStr.isEmpty())
                            bucketIdsStr += ", ";

                        bucketIdsStr += bucket;
                    }

                    if (!bucketGroup.getExtraBucketIds().isEmpty())
                    {
                        bucketIdsStr += " extra=" + bucketGroup.getExtraBucketIds().toString();
                    }
                }
                else
                {
                    // log top 10
                    List<Integer> descBucketRatioIndices = getSortedVectorIndices(bucketGroup.getBucketRatios(), false);
                    for(int j = 0; j < min(10, descBucketRatioIndices.size()); ++j)
                    {
                        Integer bucket = descBucketRatioIndices.get(j);

                        if(!bucketIdsStr.isEmpty())
                            bucketIdsStr += ", ";

                        bucketIdsStr += bucket;
                    }
                }

                final BucketFamily family = findBucketFamily(bucketGroup);

                double groupPerc = bucketGroup.getTotalCount() / totalCount;

                String linkData = "";
                if(mBucketGroupsLinked)
                {
                    linkData = String.format("family(%d) closestBg(%d css=%.4f) ",
                            family != null ? family.getId() : -1,
                            bucketGroup.getClosestBG() != null ? bucketGroup.getClosestBG().getId() : -1, bucketGroup.getClosestBGCss());
                }

                LOGGER.debug(String.format("rank %d: %s bg(%d) %scancer(%s) score(%.0f purity=%.2f LF=%.1f) samples(%d) variants(avg=%s total=%s perc=%.3f) buckets(%d: %s) effects(%s)",
                        i, groupType, bucketGroup.getId(), linkData, bucketGroup.getCancerType(), bucketGroup.calcScore(), bucketGroup.getPurity(),
                        bucketGroup.getLoadFactor(), bucketGroup.getSampleIds().size(),
                        sizeToStr(bucketGroup.getAvgCount()), sizeToStr(bucketGroup.getTotalCount()), groupPerc,
                        bucketGroup.getBucketIds().size(), bucketIdsStr, bucketGroup.getEffects()));
            }
            else
            {
                LOGGER.debug(String.format("rank %d: %s bg(%d) score(%.0f size=%d purity=%.2f) samples(%d) buckets(%d)",
                        i, groupType, bucketGroup.getId(), bucketGroup.calcScore(), bucketGroup.getSize(), bucketGroup.getPurity(),
                        bucketGroup.getSampleIds().size(), bucketGroup.getBucketIds().size()));
            }
        }
    }

    private BucketFamily findBucketFamily(final BucketGroup group)
    {
        for(final BucketFamily family : mBucketFamilies)
        {
            if(family.hasBucketGroup(group))
                return family;
        }

        return null;
    }

    private void writeBucketGroups()
    {
        try
        {
            BufferedWriter writer = getNewFile(mOutputDir, mOutputFileId + "_ba_group_data.csv");

            writer.write("Rank,BgId,FamilyId,CancerType,Effects,SampleCount,BucketCount,MutLoad,Purity");

            for(int i = 0; i < mBucketCount; ++i)
            {
                writer.write(String.format(",%d", i));
            }

            writer.newLine();

            // int minScore = 20;

            for (int bgIndex = 0; bgIndex < mBucketGroups.size(); ++bgIndex)
            {
                BucketGroup bucketGroup = mBucketGroups.get(bgIndex);

//                if (bucketGroup.getSize() < minScore)
//                    break;

                final BucketFamily family = findBucketFamily(bucketGroup);

                writer.write(String.format("%d,%d,%d,%s,%s,%d,%d,%.0f,%.2f",
                        bgIndex, bucketGroup.getId(), family != null ? family.getId() : -1, bucketGroup.getCancerType(), bucketGroup.getEffects(),
                        bucketGroup.getSampleIds().size(), bucketGroup.getBucketIds().size(),
                        sumVector(bucketGroup.getBucketCounts()), bucketGroup.getPurity()));

                double[] bucketRatios = bucketGroup.getBucketRatios();

                for(int i = 0; i < mBucketCount; ++i)
                {
                    writer.write(String.format(",%.6f", bucketRatios[i]));
                }

                writer.newLine();
            }

            writer.close();
        }
        catch(IOException exception)
        {
            LOGGER.error("failed to write output file: bucket groups");
        }
    }

    private void logSampleResults()
    {
        final List<String> categories = mExtSampleData.getFieldNames();

        int fullyAllocated = 0;
        int partiallyAllocated = 0;
        int noMatchCount = 0;
        int someMatchCount = 0;
        double countAllocated = 0; // the actual number of variants
        final double[][] elevData = mElevatedCounts.getData();

        for(int sampleId = 0; sampleId < mSampleCount; ++sampleId)
        {
            final SampleData sample = mSampleData.get(sampleId);

            final List<Integer> samBucketList = sample.getElevatedBuckets();

            if(samBucketList.isEmpty())
                continue;

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

            if(sample.getAllocPercent() >= MIN_BUCKET_ALLOCATION)
            {
                ++fullyAllocated;
                LOGGER.debug(String.format("sample(%d: %s) fully allocated: groups(%d) buckets(%d) count(%s alloc=%.3f) cancer(%s) effects(%d: %s)",
                        sampleId, sample.getSampleName(), sample.getElevBucketGroups().size(), samBucketList.size(),
                        sizeToStr(sample.getElevatedCount()), sample.getAllocPercent(), cancerType, effectsCount, effects));
                continue;
            }

            if(sample.getAllocPercent() > 0.1)
            {
                ++partiallyAllocated;
                LOGGER.debug(String.format("sample(%d: %s) partially allocated: groups(%d) buckets(%d) count(%s alloc=%.3f) cancer(%s) effects(%d: %s)",
                        sampleId, sample.getSampleName(), sample.getElevBucketGroups().size(), samBucketList.size(),
                        sizeToStr(sample.getElevatedCount()), sample.getAllocPercent(), cancerType, effectsCount, effects));
                continue;
            }

            // now find the closest matching bucket(s) and try to work out why they didn't match
            boolean someMatch = false;
            for(final BucketGroup bucketGroup : mBucketGroups)
            {
                if(!bucketGroup.getCancerType().equals(cancerType))
                    continue;

                final List<Integer> commonBuckets = getMatchingList(bucketGroup.getBucketIds(), samBucketList);
                int commonBucketCount = commonBuckets.size();

                if(commonBucketCount < 0.5 * samBucketList.size() || commonBucketCount < 2)
                    continue;

                double[] sc1 = extractBucketCountSubset(sampleId, commonBuckets);
                double bcCss = calcSharedCSS(sc1, bucketGroup.getBucketCounts());

                if(bcCss >= 0.9)
                {
                    someMatch = true;

                    LOGGER.debug(String.format("sample(%d: %s ct=%s eff=%d: %s) unallocated close match with bg(%d eff=%s), buckets(grp=%d sam=%d match=%d) css(%.4f)",
                            sampleId, sample.getSampleName(), cancerType, effectsCount, effects, bucketGroup.getId(), bucketGroup.getEffects(),
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
                        sampleId, sample.getSampleName(), samBucketList.size(), cancerType, effectsCount, effects);
                ++noMatchCount;
            }
        }

        int unallocated = someMatchCount + noMatchCount;
        double percAllocated = countAllocated / mElevatedCount;

        LOGGER.debug(String.format("sample summary: total(%d) elev(%d) alloc(%d, %.3f of %s) partial(%d) unalloc(%d possibleMatch=%d noMatch=%d) groupCount(%d)",
                mSampleCount, mAllSampleBucketGroups.size(), fullyAllocated, percAllocated, sizeToStr(mElevatedCount),
                partiallyAllocated, unallocated, someMatchCount, noMatchCount, mBucketGroups.size()));
    }

    private void compareSignatures()
    {
        double sigCompareCss = 0.95;

        if(mProposedSigs != null)
        {
            // first the internally generated ones
            List<double[]> cssResults = getTopCssPairs(mProposedSigs, mProposedSigs, sigCompareCss, true, true);

            if (cssResults.isEmpty())
            {
                LOGGER.debug("no similar proposed sigs from bucket groups");
            } else
            {
                for (final double[] result : cssResults)
                {
                    int sigId1 = (int)result[CSSR_I1];
                    int sigId2 = (int)result[CSSR_I2];
                    final BucketGroup bg1 = mBucketGroups.get(mSigToBgMapping.get(sigId1));
                    final BucketGroup bg2 = mBucketGroups.get(mSigToBgMapping.get(sigId2));

                    LOGGER.debug(String.format("proposed sig(%s bg=%d: ct=%s eff=%s samples=%d) matches sig(%d bg=%d: ct=%s eff=%s samples=%d) with css(%.4f)",
                            sigId1, bg1.getId(), bg1.getCancerType(), bg1.getEffects(), bg1.getSampleIds().size(),
                            sigId2, bg2.getId(), bg2.getCancerType(), bg2.getEffects(), bg2.getSampleIds().size(), result[CSSR_VAL]));
                }
            }
        }

        if(mReferenceSigs == null)
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
                int externalSigId = (int)result[CSSR_I1] + 1; // bumped up to correspond to convention of starting with 1
                int proposedSigId = (int)result[CSSR_I2];
                final BucketGroup bucketGroup = mBucketGroups.get(mSigToBgMapping.get(proposedSigId));

                LOGGER.debug(String.format("external ref sig(%d) matches bg(%d: ct=%s eff=%s samples=%d purity=%.2f) with css(%.4f)",
                        externalSigId, bucketGroup.getId(), bucketGroup.getCancerType(), bucketGroup.getEffects(),
                        bucketGroup.getSampleIds().size(), bucketGroup.getPurity(), result[CSSR_VAL]));
            }
        }
    }

    private void formBucketFamilies(boolean useElevated)
    {
        List<BucketGroup> bgList = useElevated ? mBucketGroups : mBackgroundGroups;

        Collections.sort(bgList);

        LOGGER.debug("creating families from {} {} bucket groups", bgList.size(), useElevated ? "elevated" : "background");

        int groupCount = min(bgList.size(), 3000);

        NmfMatrix bgMatrix = new NmfMatrix(mBucketCount, groupCount);

        for(int bgIndex = 0; bgIndex < groupCount; ++bgIndex)
        {
            final BucketGroup bg = bgList.get(bgIndex);
            bgMatrix.setCol(bgIndex, bg.getBucketRatios());
        }

        bgMatrix.cacheTranspose();

        // use a lower CSS check to pick up similiarities to store against the closest-match data on the BG
        List<double[]> cssResults = getTopCssPairs(bgMatrix, bgMatrix, 0.95, false, true);

        if (cssResults.isEmpty())
        {
            LOGGER.debug("no similar bucket groups to create families");
            return;
        }

        List<Integer> assignedBGs = Lists.newArrayList();

        for (int r1 = 0; r1 < cssResults.size(); ++r1)
        {
            final double[] result1 = cssResults.get(r1);

            int bg1 = (int) result1[CSSR_I1];
            int bg2 = (int) result1[CSSR_I2];
            double css = result1[CSSR_VAL];

            // keep track of closest match for analysis purposes
            final BucketGroup bucketGroup1 = bgList.get(bg1);
            final BucketGroup bucketGroup2 = bgList.get(bg2);

            if(bucketGroup1.getClosestBG() == null)
            {
                bucketGroup1.setClosestBG(bucketGroup2);
                bucketGroup1.setClosestBGCss(css);
            }

            if(bucketGroup2.getClosestBG() == null)
            {
                bucketGroup2.setClosestBG(bucketGroup1);
                bucketGroup2.setClosestBGCss(css);
            }

            if(css < mSigCssThreshold)
                continue;

            if (assignedBGs.contains(bg1) || assignedBGs.contains(bg2))
                continue;

            List<Integer> similarBGs = Lists.newArrayList();

            double cssTotal = result1[CSSR_VAL];
            int cssMatchCount = 1;

            similarBGs.add(bg1);
            similarBGs.add(bg2);

            for (int r2 = r1 + 1; r2 < cssResults.size(); ++r2)
            {
                final double[] result2 = cssResults.get(r2);

                int bg3 = (int) result2[CSSR_I1];
                int bg4 = (int) result2[CSSR_I2];

                if (assignedBGs.contains(bg3) || assignedBGs.contains(bg3))
                    continue;

                if (bg1 == bg3 || bg1 == bg4 || bg2 == bg3 || bg2 == bg4)
                {
                    if (!similarBGs.contains(bg3))
                        similarBGs.add(bg3);

                    if (!similarBGs.contains(bg4))
                        similarBGs.add(bg4);

                    cssTotal += result1[CSSR_VAL];
                    ++cssMatchCount;
                }
            }

            BucketFamily bucketFamily = new BucketFamily(mBucketFamilies.size());

            for(Integer bgIndex : similarBGs)
            {
                final BucketGroup bucketGroup = bgList.get(bgIndex);

                bucketFamily.addBucketGroup(bucketGroup);

                // keep track of which groups have been added to a family
                assignedBGs.add(bgIndex);
            }

            mBucketFamilies.add(bucketFamily);

            bucketFamily.calcAll(mElevatedCounts);

            double totalCount = useElevated ? mElevatedCount : mBackgroundCounts.sum();
            double familyPerc = bucketFamily.getTotalCount() / totalCount;

            // final BucketGroup bucketGroup = mBucketGroups.get(mSigToBgMapping.get(proposedSigId));
            LOGGER.debug(String.format("%s family(%d) created from %d bucketGroups, buckets(%d) samples(%d) count(%s perc=%.3f) avgCss(%.4f) cancer(%s) effects(%s)",
                    useElevated ? "elevated" : "background", bucketFamily.getId(), bucketFamily.getBucketGroups().size(), bucketFamily.getBucketIds().size(), bucketFamily.getSampleIds().size(),
                    sizeToStr(bucketFamily.getTotalCount()), familyPerc, cssTotal/cssMatchCount, bucketFamily.getCancerType(), bucketFamily.getEffects()));
        }

        mBucketGroupsLinked = true;
    }

    private void formBackgroundBucketGroups()
    {
        for (Map.Entry<String, List<Integer>> entry : mCancerSamplesMap.entrySet())
        {
            final String cancerType = entry.getKey();
            List<Integer> sampleIds = entry.getValue();

            assignToBackgroundBucketGroups(cancerType, sampleIds);
        }
    }

    private void assignToBackgroundBucketGroups(final String cancerType, final List<Integer> sampleIds)
    {
        LOGGER.debug("cancerType({}) creating background groups for {} samples", cancerType, sampleIds.size());

        List<Integer> fullBucketSet = Lists.newArrayList();
        for (int i = 0; i < mBucketCount; ++i)
        {
            fullBucketSet.add(i);
        }

        BucketGroup bucketGroup = new BucketGroup(mBackgroundGroups.size());
        bucketGroup.addBuckets(fullBucketSet);
        bucketGroup.setCancerType(cancerType);
        mBackgroundGroups.add(bucketGroup);

        final List<Double> medianCounts = mBucketMediansMap.get(cancerType);
        final double[][] bgData = mBackgroundCounts.getData();

        if(medianCounts == null)
            return;

        double[] bucketRatios = listToArray(medianCounts);

        // for now add all samples to a single group per cancer type
        for (int index = 0; index < sampleIds.size(); ++index)
        {
            int sampleId = sampleIds.get(index);
            double[] sampleCounts = mBackgroundCounts.getCol(sampleId);

            /*
            double bcCss = calcCSS(sampleCounts, bucketRatios);

            boolean withinProb = isWithinProbability(sampleId, bucketRatios, sampleCounts, false);
            double sampleTotal = sumVector(sampleCounts);

            int lowProbCount = 0;
            int highProbCount = 0;

            if (!withinProb)
            {
                int[] result = checkPermittedSampleCounts(sampleId, bucketRatios, sampleCounts, false, false);
                lowProbCount = result[PI_LOW];
                highProbCount = result[PI_HIGH];
            }
            */

            // add regardless of the match
            bucketGroup.addSample(sampleId, sampleCounts);

//            LOGGER.debug(String.format("bg(%d) added sample(%d) varCount(%.0f) css(%.4f prob=%s low=%d high=%d)",
//                    bucketGroup.getId(), sampleId, sampleTotal, bcCss, withinProb ? "pass" : "fail",
//                    lowProbCount, highProbCount));
        }

        // now counts are set, manually set the ratios (ie irrespective of the sample counts)
        bucketGroup.setBucketRatios(bucketRatios);
    }

    private void formBackgroundBucketGroups(final String cancerType, final List<Integer> sampleIds)
    {
        LOGGER.debug("cancerType({}) creating background groups for {} samples", cancerType, sampleIds.size());

        int samplesAdded = 0;

        List<BucketGroup> bgList = Lists.newArrayList();

        List<Integer> fullBucketSet = Lists.newArrayList();
        for(int i = 0; i < mBucketCount; ++i)
        {
            fullBucketSet.add(i);
        }

        for (int index1 = 0; index1 < sampleIds.size(); ++index1)
        {
            int samIndex1 = sampleIds.get(index1);
            double[] sc1 = mBackgroundCounts.getCol(samIndex1);

            boolean allocated = false;

            for (BucketGroup bucketGroup : bgList)
            {
                if(bucketGroup.hasSample(samIndex1))
                {
                    allocated = true;
                    break;
                }

                double bcCss = calcCSS(sc1, bucketGroup.getBucketCounts());

                boolean withinProb = isWithinProbability(samIndex1, bucketGroup.getBucketRatios(), sc1, false);
                boolean cssOK = bcCss >= mHighCssThreshold;

                if(cssOK || withinProb)
                {
                    if(cssOK && !withinProb)
                    {
                        int[] result = checkPermittedSampleCounts(samIndex1, bucketGroup.getBucketRatios(), sc1, false, false);

                        double lowPerc = result[PI_LOW] / (double)result[PI_COUNT];
                        if(lowPerc >= 0.05)
                        {
                            double sampleBucketCount = sumVector(sc1);
                            LOGGER.debug(String.format("bg(%d) sample(%d) varCount(%.0f) has css(%.4f) prob(low=%d high=%d) mismatch",
                                    bucketGroup.getId(), samIndex1, sampleBucketCount, bcCss, result[PI_LOW], result[PI_HIGH]));
                            continue;
                        }
                    }

                    // group's buckets remain the same, neither increased nor refined, just add this new sample's counts
                    bucketGroup.addSample(samIndex1, sc1);
                    ++samplesAdded;

                    LOGGER.debug(String.format("bg(%d) added sample(%d) css(%.4f prob=%s) totalSamples(%d)",
                            bucketGroup.getId(), samIndex1, bcCss, withinProb ? "pass" : "fail", bucketGroup.getSampleIds().size()));

                    allocated = true;
                    break;
                }
            }

            if (allocated)
            {
                continue;
            }

            for (int index2 = index1 + 1; index2 < sampleIds.size(); ++index2)
            {
                int samIndex2 = sampleIds.get(index2);
                double[] sc2 = mBackgroundCounts.getCol(samIndex2);
                double bcCss = calcCSS(sc1, sc2);

                if(bcCss >= mHighCssThreshold)
                {
                    BucketGroup bucketGroup = new BucketGroup(mBackgroundGroups.size() + bgList.size());
                    bucketGroup.addBuckets(fullBucketSet);
                    bucketGroup.addSample(samIndex1, sc1);
                    bucketGroup.addSample(samIndex2, sc2);
                    samplesAdded += 2;

                    LOGGER.debug(String.format("added bg(%d) samples(%d and %d) css(%.4f)",
                            bucketGroup.getId(), samIndex1, samIndex2, bcCss));

                    bgList.add(bucketGroup);
                    break;
                }
            }
        }

        if(bgList.isEmpty())
        {
            LOGGER.debug("no background groups created");
            return;
        }

        mBackgroundGroups.addAll(bgList);

        LOGGER.debug("cancerType({}) created {} background groups with {} samples, total groups({})",
                cancerType, bgList.size(), samplesAdded, mBackgroundGroups.size());
    }

    private void writeBucketFamilies()
    {
        try
        {
            BufferedWriter writer = getNewFile(mOutputDir, mOutputFileId + "_ba_families.csv");

            writer.write("FamilyId,BgId,CancerType,Effects,SampleCount,BucketCount,MutLoad,Purity");

            for (int i = 0; i < mBucketCount; ++i)
            {
                writer.write(String.format(",%d", i));
            }

            writer.newLine();

            for(BucketFamily bucketFamily : mBucketFamilies)
            {
                for(BucketGroup bucketGroup : bucketFamily.getBucketGroups())
                {
                    writer.write(String.format("%d,%d,%s,%s,%d,%d,%.0f,%.2f",
                            bucketFamily.getId(), bucketGroup.getId(), bucketGroup.getCancerType(), bucketGroup.getEffects(),
                            bucketGroup.getSampleIds().size(), bucketGroup.getBucketIds().size(),
                            sumVector(bucketGroup.getBucketCounts()), bucketGroup.getPurity()));

                    double[] bucketRatios = bucketGroup.getBucketRatios();

                    for(int i = 0; i < mBucketCount; ++i)
                    {
                        writer.write(String.format(",%.6f", bucketRatios[i]));
                    }

                    writer.newLine();
                }
            }

            writer.close();
        }
        catch (IOException exception)
        {
            LOGGER.error("failed to write output file: bucket families");
        }
    }

    private void writeBucketGroupOverlap()
    {
        try
        {
            BufferedWriter writer = getNewFile(mOutputDir, mOutputFileId + "_ba_group_overlap.csv");

            writer.write("BgId1,SC1,BC1,BgId2,SC2,BC2,SharedBuckets,SharedSamples,OverlapSize,Perc1Of2,Perc2Of1");
            writer.newLine();

            for(int bgIndex1 = 0; bgIndex1 < mBucketGroups.size(); ++bgIndex1)
            {
                final BucketGroup bg1 = mBucketGroups.get(bgIndex1);

                for (int bgIndex2 = bgIndex1 + 1; bgIndex2 < mBucketGroups.size(); ++bgIndex2)
                {
                    final BucketGroup bg2 = mBucketGroups.get(bgIndex2);

                    final List<Integer> bl1 = bg1.getBucketIds();
                    final List<Integer> bl2 = bg2.getBucketIds();

                    final List<Integer> sharedBuckets = getMatchingList(bl1, bl2);

                    if(sharedBuckets.isEmpty())
                        continue;

                    final List<Integer> sl1 = bg1.getSampleIds();
                    final List<Integer> sl2 = bg2.getSampleIds();

                    final List<Integer> sharedSamples = getMatchingList(sl1, sl2);

                    int sharedSize = sharedBuckets.size() * sharedSamples.size();

                    if(sharedSize == 0)
                        continue;

                    writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%d,%d,%.3f,%.3f",
                            bg1.getId(), sl1.size(), bl1.size(), bg2.getId(),
                            sl2.size(), bl2.size(), sharedBuckets.size(), sharedSamples.size(),
                            sharedSize, sharedSize / (double)bg1.getSize(), sharedSize / (double)bg2.getSize()));

                    writer.newLine();
                }
            }

            writer.close();
        }
        catch(IOException exception)
        {
            LOGGER.error("failed to write output file: bucket groups");
        }
    }

    private void createSignatures()
    {
        LOGGER.debug("creating signatures");

        int proposedSigCount = min(mMaxProposedSigs, mBucketGroups.size());

        mProposedSigs = new NmfMatrix(mBucketCount, proposedSigCount);

        NmfMatrix sigTest = new NmfMatrix(mBucketCount, 1);

        double purityThreshold = 0.7;
        int scoreThreshold = 50; // eg 100% pure, 2 samples 10 buckets, LF=2.5
        int minSamples = 2;

        int sigId = 0;

        for(int bgIndex = 0; bgIndex < mBucketGroups.size(); ++bgIndex)
        {
            final BucketGroup bucketGroup = mBucketGroups.get(bgIndex);

            if(bucketGroup.getPurity() < purityThreshold || bucketGroup.calcScore() < scoreThreshold
            || bucketGroup.getSampleIds().size() < minSamples)
                continue;

            final double[] bucketRatios = bucketGroup.getBucketRatios();

            sigTest.setCol(0, bucketRatios);

            List<double[]> cssResults = getTopCssPairs(sigTest, mProposedSigs, mSigCssThreshold, false, false);

            if (!cssResults.isEmpty())
            {
                final double[] result = cssResults.get(0);
                // LOGGER.debug(String.format("bg(%d) too similar to sig(%d) with css(%.4f)", bucketGroup.getId(), (int)result[CSSR_I2], result[CSSR_VAL]));
                continue;
            }

            LOGGER.debug(String.format("bg(%d) added proposed sig(%d): score(%.0f purity=%.2f sam=%d buckets=%s) type(%s effects=%s)",
                    bucketGroup.getId(), sigId, bucketGroup.calcScore(), bucketGroup.getPurity(), bucketGroup.getSampleIds().size(),
                    bucketGroup.getBucketIds().size(), bucketGroup.getCancerType(), bucketGroup.getEffects()));

            mProposedSigs.setCol(sigId, bucketRatios);
            mSigToBgMapping.add(bgIndex);

            ++sigId;

            if(sigId >= proposedSigCount)
                break;
        }

        if(sigId != proposedSigCount)
        {
            mProposedSigs = redimension(mProposedSigs, mProposedSigs.Rows, sigId);
        }

        writeSignatures(mProposedSigs, "_proposed_sigs.csv", null);
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

    private List<Integer> getSampleBucketGroupIds(int sampleId)
    {
        return getSampleBucketGroups(sampleId, true);
    }

    private List<Integer> getSampleBucketGroupIndices(int sampleId)
    {
        return getSampleBucketGroups(sampleId, false);
    }

    private List<Integer> getSampleBucketGroups(int sampleId, boolean getIds)
    {
        List<Integer> groupIds = Lists.newArrayList();

        for (int i = 0; i < mBucketGroups.size(); ++i)
        {
            BucketGroup bucketGroup = mBucketGroups.get(i);

            if (bucketGroup.hasSample(sampleId))
                groupIds.add(getIds ? bucketGroup.getId() : i);
        }

        return groupIds;
    }

    private int getSampleGroupMembershipCount(int sample)
    {
        return getSampleBucketGroupIndices(sample).size();
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
        double[] vec = new double[mBucketCount];

        // only the values above the expected count are used for sample-count comparisons
        final double[][] elevData = mElevatedCounts.getData();

        for(int i = 0; i < bucketSubset.size(); ++i)
        {
            int bucketIndex = bucketSubset.get(i);
            vec[bucketIndex] = elevData[bucketIndex][sam1];
        }

        return vec;
    }

    private double calcSharedCSS(final double[] set1, final double[] set2)
    {
        return calcCSS(set1, set2, true);
        // return calcCSSRelative(set1, set2);
    }

    private boolean isWithinProbability(int sampleId, final double[] ratios, final double[] sampleData, boolean useElevated)
    {
        int[] result = checkPermittedSampleCounts(sampleId, ratios, sampleData, useElevated, false);
        return result[PI_LOW] == 0 && result[PI_HIGH] == 0;
    }

    private static int PI_COUNT = 0;
    private static int PI_LOW = 1;
    private static int PI_HIGH = 2;

    private int[] checkPermittedSampleCounts(int sampleId, final double[] ratios, final double[] sampleData, boolean useElevated, boolean verbose)
    {
        double sampleTotal = sumVector(sampleData);

        final double[][] prData = useElevated ? mPermittedElevRange.getData() : mPermittedBgRange.getData();
        final double[][] countData = useElevated ? mElevatedCounts.getData() : mBackgroundCounts.getData();

        int[] result = new int[3];

        for(int i = 0; i < sampleData.length; ++i)
        {
            if(sampleData[i] == 0)
                continue;

            if(prData[i][sampleId] == 0)
                continue;

            ++result[PI_COUNT];

            double bucketValue = ratios[i] * sampleTotal;
            double rangeHigh = countData[i][sampleId] + prData[i][sampleId];
            double rangeLow = countData[i][sampleId] - prData[i][sampleId];

            if(verbose && (bucketValue > rangeHigh || bucketValue < rangeLow))
            {
                LOGGER.debug("sample({}) value({}) vs range({} -> {})", sampleId, round(bucketValue), rangeLow, rangeHigh);
            }

            if(bucketValue > rangeHigh)
                ++result[PI_HIGH];
            else if(bucketValue < rangeLow)
                ++result[PI_LOW];
        }

        // double llProb = calcLogLikelihood(sampleData, impliedCounts, false);

        return result;
    }

    private void populateCancerSamplesMap()
    {
        if(mSampleData.isEmpty())
            return;

        mCancerSamplesMap = new HashMap();

        String minors = "Minors";
        int minSamples = 20;

        List<String> cancerTypes = Lists.newArrayList();

        for(final SampleData sample : mSampleData)
        {
            final String cancerType = sample.getCancerType();
            List<Integer> samplesList = mCancerSamplesMap.get(cancerType);

            if(samplesList == null)
            {
                samplesList = Lists.newArrayList();
                mCancerSamplesMap.put(cancerType, samplesList);
                cancerTypes.add(cancerType);
            }

            samplesList.add(sample.Id);
        }


        List<Integer> minorsList = Lists.newArrayList();
        mCancerSamplesMap.put(minors, minorsList);

        for(final String cancerType : cancerTypes)
        {
            List<Integer> samplesList = mCancerSamplesMap.get(cancerType);

            if(samplesList.size() < minSamples)
            {
                minorsList.addAll(samplesList);
                mCancerSamplesMap.remove(cancerType);
            }
        }
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

    private void writeBackgroundSigs()
    {
        int sigCount = mBackgroundGroups.size();

        NmfMatrix bgSigs = new NmfMatrix(mBucketCount, sigCount);

        List<String> cancerTypes = Lists.newArrayList();

        for(int sig = 0; sig < sigCount; ++sig)
        {
            final BucketGroup bucketGroup = mBackgroundGroups.get(sig);
            bgSigs.setCol(sig, bucketGroup.getBucketRatios());
            cancerTypes.add(bucketGroup.getCancerType());
        }

        writeSignatures(bgSigs, "_ba_background_sigs.csv", cancerTypes);
    }

    private void writeSampleData()
    {
        try
        {
            BufferedWriter writer = getNewFile(mOutputDir,mOutputFileId + "_ba_sample_alloc.csv");

            writer.write("SampleIndex,SampleId,MutLoad,BackgroundCount,ElevCount,Status,AllocPerc,BgIds,BgCancerType,BgEffects,Buckets");
            writer.newLine();

            final double[][] elevData = mElevatedCounts.getData();
            final double[][] backgroundData = mBackgroundCounts.getData();

            for(int i = 0; i < mSampleCount; ++i)
            {
                final SampleData sampleData = mSampleData.get(i);

                final List<Integer> sampleBuckets = sampleData.getElevatedBuckets();

                writer.write(String.format("%d,%s", i, sampleData.getSampleName()));

                double bgTotal = sumVector(mBackgroundCounts.getCol(i));

                writer.write(String.format(",%.0f,%.0f,%.0f", mSampleTotals[i], bgTotal, sampleData.getElevatedCount()));

                if(sampleBuckets == null)
                {
                    writer.write(",NoElevation,,,,,");
                    writer.newLine();
                    continue;
                }

                double allocPerc = sampleData.getAllocPercent();

                String status = "Unalloc";
                if(allocPerc >= 0.9)
                {
                    status = "Alloc";
                }
                else if(allocPerc > 0.1)
                {
                    status = "Partial";
                }

                writer.write(String.format(",%s,%.3f", status, allocPerc));

                final List<BucketGroup> bucketGroups = sampleData.getElevBucketGroups();
                if(!bucketGroups.isEmpty())
                {
                    String bgIdsStr = "";

                    for(final BucketGroup bucketGroup : bucketGroups)
                    {
                        if(!bgIdsStr.isEmpty())
                            bgIdsStr += ";";

                        bgIdsStr += bucketGroup.getId();
                    }

                    writer.write(String.format(",%s,%s,%s",
                            bgIdsStr, bucketGroups.get(0).getCancerType(), bucketGroups.get(0).getEffects()));
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

    public void writeSignatures(final NmfMatrix signatures, final String fileId, final List<String> sigIds)
    {
        try
        {
            BufferedWriter writer = getNewFile(mOutputDir,mOutputFileId + fileId); // eg "_ba_sigs.csv");

            if(sigIds != null)
            {
                for(int i = 0; i < sigIds.size(); ++i)
                {
                    if(i > 0)
                        writer.write(String.format(",%s", sigIds.get(i)));
                    else
                        writer.write(String.format("%s", sigIds.get(i)));
                }
            }
            else
            {
                int i = 0;
                for (; i < signatures.Cols - 1; ++i)
                {
                    writer.write(String.format("%d,", i));
                }
                writer.write(String.format("%d", i));
            }

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


    public void writeSampleCalcData()
    {
        if(mLoadedSampleCalcData)
            return;

        try
        {
            BufferedWriter writer = getNewFile(mOutputDir, mOutputFileId + "_ba_sam_calc_data.csv");

            writer.write("SampleId,SampleName,BgAllocation");
            writer.newLine();

            for(int sampleId = 0; sampleId < mSampleCount; ++sampleId)
            {
                final SampleData sample = mSampleData.get(sampleId);

                writer.write(String.format("%d,%s", sampleId, sample.getSampleName()));

                // computed fields to cache
                writer.write(String.format(",%.0f", mSampleBgAllocations.get(sampleId)));

                writer.newLine();;
            }

            writer.close();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile");
        }
    }


            // old methods
    /*
    private void formSampleSubGroups()
    {
        // forms groups out of samples with similarly elevated buckets
        // returns true if groups were created or added to
        int groupsAdjusted = 0;
        int groupsCreated = 0;

        for (int samIndex1 = 0; samIndex1 < mSampleCount; ++samIndex1)
        {
            final List<Integer> bl1 = mWorkingSBGroups.get(samIndex1);

            if (bl1 == null)
                continue;

            if(mSampleWatchList.contains(samIndex1))
            {
                LOGGER.debug("specific sample");
            }

            boolean sam1Allocated = false;

            for (int samIndex2 = samIndex1 + 1; samIndex2 < mSampleCount; ++samIndex2)
            {
                final List<Integer> bl2 = mWorkingSBGroups.get(samIndex2);

                if (bl2 == null)
                    continue;

                final List<Integer> commonBuckets = getMatchingList(bl1, bl2);
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

        int groupsAdjusted = 0;
        int groupsCreated = 0;

        for (int samIndex1 = 0; samIndex1 < mSampleCount; ++samIndex1)
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
                final List<Integer> commonBuckets = getMatchingList(groupBuckets, bl1);

                if(areAnySampleBucketsAllocated(samIndex1, commonBuckets))
                {
                    continue;
                }

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

                    final List<Integer> commonBuckets = getMatchingList(bucketGroup.getBucketIds(), bl1);

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

            for (int samIndex2 = samIndex1 + 1; samIndex2 < mSampleCount; ++samIndex2)
            {
                final List<Integer> bl2 = mWorkingSBGroups.get(samIndex2);

                if (bl2 == null)
                    continue;

                final List<Integer> commonBuckets = getMatchingList(bl1, bl2);

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
    */

    private void calcSampleOverlaps()
    {
        LOGGER.debug("calculating sample overlaps");

        try
        {
            BufferedWriter writer = getNewFile(mOutputDir,mOutputFileId + "_ba_all_sample_css.csv");

            writer.write("SamId1,SamName2,SamCT1,S1Elevated,SamId2,SamName2,SamCT2,S2Elevated,Common,Matched,CSS,Buckets");
            writer.newLine();

            final double[][] scData = mElevatedCounts.getData();

            for (int samIndex1 = 0; samIndex1 < mSampleCount; ++samIndex1)
            {
                final List<Integer> bl1 = mAllSampleBucketGroups.get(samIndex1);

                if (bl1 == null)
                    continue;

                final SampleData sample1 = mSampleData.get(samIndex1);

                for (int samIndex2 = samIndex1 + 1; samIndex2 < mSampleCount; ++samIndex2)
                {
                    final List<Integer> bl2 = mAllSampleBucketGroups.get(samIndex2);

                    if (bl2 == null)
                        continue;

                    List<Integer> commonBuckets = getMatchingList(bl1, bl2);
                    int commonBucketCount = commonBuckets.size();

                    if(commonBucketCount < 2)
                        continue;

                    double[] sc1 = extractBucketCountSubset(samIndex1, commonBuckets);
                    double[] sc2 = extractBucketCountSubset(samIndex2, commonBuckets);
                    double bcCss = calcSharedCSS(sc1, sc2);

                    double cssMatched = 0;
                    int cssMatchedCount = 0;

                    List<Integer> removedBuckets = Lists.newArrayList();

                    if(bcCss >= mHighCssThreshold)
                    {
                        cssMatched = bcCss;
                        cssMatchedCount = commonBucketCount;
                    }
                    else if(commonBucketCount > 2)
                    {
                        // attempt to find a match using less overlapping buckets
                        int currentBucketCount = commonBucketCount;

                        while(currentBucketCount > 2 && bcCss < mHighCssThreshold)
                        {
                            int lastTestBucket = -1;
                            double bestCss = 0;
                            int bestTestBucket = 0;

                            for(Integer removedBucket : removedBuckets)
                            {
                                sc1[removedBucket] = 0;
                                sc2[removedBucket] = 0;
                            }

                            for(Integer testBucket : commonBuckets)
                            {
                                if(removedBuckets.contains(testBucket)) // no point in trying these each time
                                    continue;

                                sc1[testBucket] = 0;
                                sc2[testBucket] = 0;

                                // restore previous test bucket
                                if(lastTestBucket != -1)
                                {
                                    sc1[lastTestBucket] = scData[lastTestBucket][samIndex1];
                                    sc2[lastTestBucket] = scData[lastTestBucket][samIndex2];
                                }

                                // run CSS on this reduced set of buckets
                                bcCss = calcSharedCSS(sc1, sc2);
                                lastTestBucket = testBucket;

                                if(bcCss > bestCss)
                                {
                                    bestCss = bcCss;
                                    bestTestBucket = testBucket;
                                }

                                if(bcCss >= mHighCssThreshold)
                                    break;
                            }

                            // restore last test bucket
                            if(lastTestBucket != -1)
                            {
                                sc1[lastTestBucket] = scData[lastTestBucket][samIndex1];
                                sc2[lastTestBucket] = scData[lastTestBucket][samIndex2];
                            }

                            --currentBucketCount;
                            cssMatched = bestCss;
                            removedBuckets.add(bestTestBucket);

                            if(bestCss >= mHighCssThreshold)
                            {
                                cssMatchedCount = currentBucketCount;
                            }
                        }
                    }

                    String bucketIdsStr = "";
                    for(Integer bucket : commonBuckets)
                    {
                        if(removedBuckets.contains(bucket))
                            continue;

                        if(!bucketIdsStr.isEmpty())
                            bucketIdsStr += ";";

                        bucketIdsStr += bucket;
                    }

                    final SampleData sample2 = mSampleData.get(samIndex2);

                    writer.write(String.format("%d,%s,%s,%d,%d,%s,%s,%d,%d,%d,%.6f,%s",
                            samIndex1, sample1.getSampleName(), sample1.getCancerType(), bl1.size(),
                            samIndex2, sample2.getSampleName(), sample2.getCancerType(), bl2.size(),
                            commonBucketCount, cssMatchedCount, cssMatched, bucketIdsStr));

                    writer.newLine();
                }
            }

            writer.close();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile");
        }
    }


}
