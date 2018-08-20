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
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.calcMinPositiveRatio;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.capValue;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.convertToPercentages;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.copyVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getCombinedList;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getDiffList;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getMatchingList;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getNewFile;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.listToArray;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sizeToStr;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sumVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.vectorMultiply;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.writeMatrixData;
import static com.hartwig.hmftools.data_analyser.calcs.NmfConfig.NMF_REF_SIG_FILE;
import static com.hartwig.hmftools.data_analyser.calcs.NmfSampleFitter.fitCountsToRatios;
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
    private double mElevatedCount;
    private double mAllocatedCount;
    private double mBackgroundCount;
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

    private int mNextBucketId;
    private List<BucketGroup> mBucketGroups;
    private List<BucketGroup> mFinalBucketGroups;
    private List<BucketGroup> mTopAllocBucketGroups;
    private List<BucketGroup> mBackgroundGroups;
    private List<BucketFamily> mBucketFamilies;
    private boolean mBucketGroupsLinked;
    private List<Integer> mSkippedSamples;
    private List<Integer> mReassessSamples;

    private boolean mNoBackgroundCounts; // whether to make a distinction between background and elevated counts
    private boolean mApplyNoise; // whether to factor Poisson noise into the sample counts and fits

    BufferedWriter mBucketGroupFileWriter;

    // config
    private String mOutputDir;
    private String mOutputFileId;
    private String mSpecificCancer;
    private double mHighCssThreshold; // CSS level for samples or groups to be consider similar
    private double mHighRequiredMatch; // high-level required number of buckets to match between samples or groups
    private int mMaxProposedSigs;
    private double mSigCssThreshold; // avoid creating sigs with similarity lower than this
    private int mMutationalLoadCap;
    private boolean mLogVerbose;

    private static String BA_CSS_HIGH_THRESHOLD = "ba_css_high";
    private static String BA_MAX_PROPOSED_SIGS = "ba_max_proposed_sigs";
    private static String BA_CSS_SIG_THRESHOLD = "ba_css_proposed_sigs";
    private static String BA_SPECIFIC_CANCER = "ba_specific_cancer";

    // constraint constants - consider moing any allocation-related ones to config to make them visible
    private static double SAMPLE_ALLOCATED_PERCENT = 0.995;
    private static double DOMINANT_CATEGORY_PERCENT = 0.7; /// mark a group with a category if X% of samples in it have this attribute (eg cancer type, UV)
    private static double MAX_ELEVATED_PROB = 1e-12;
    private int MIN_BUCKET_COUNT_OVERLAP = 3; // used for pairings of samples with reduced bucket overlap
    private static double PERMITTED_PROB_NOISE = 1e-4;
    private static double MIN_GROUP_ALLOC_PERCENT = 0.1; // only allocate a sample to a group if it takes at this much of the elevated count
    private static double MIN_GROUP_ALLOC_PERCENT_LOWER = 0.03; // hard lower limit
    private static double SKIP_ALLOC_FACTOR = 2;

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
        mAllocatedCount = 0;
        mBackgroundCount = 0;
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

        mNextBucketId = 0;
        mBucketGroups = Lists.newArrayList();
        mTopAllocBucketGroups = Lists.newArrayList();
        mFinalBucketGroups = Lists.newArrayList();
        mBackgroundGroups = Lists.newArrayList();
        mBucketFamilies = Lists.newArrayList();
        mSkippedSamples = Lists.newArrayList();
        mReassessSamples = Lists.newArrayList();
        mBucketGroupsLinked = false;
        mBucketGroupFileWriter = null;

        mHighCssThreshold = 0.995;
        mHighRequiredMatch = 0.95;
        mMaxProposedSigs = 0;
        mSigCssThreshold = mHighCssThreshold * 0.95;
        mMutationalLoadCap = 0;
        mLogVerbose = false;

        mNoBackgroundCounts = true;
        mApplyNoise = true;

        mSampleWatchList = Lists.newArrayList();

        mSampleWatchList.add(2075);
        // mSampleWatchList.add(1714);
        //mSampleWatchList.add(1375);
//        mSampleWatchList.add(2617);
    }

    public static void addCmdLineArgs(Options options)
    {

        options.addOption(BA_EXT_SAMPLE_DATA_FILE, true, "Sample external data");
        options.addOption(BA_CSS_HIGH_THRESHOLD, true, "Cosine sim for high-match test");
        options.addOption(BA_CSS_SIG_THRESHOLD, true, "Cosine sim for comparing proposed sigs");
        options.addOption(BA_MAX_PROPOSED_SIGS, true, "Maximum number of bucket groups to turn into proposed sigs");
        options.addOption(BA_SAMPLE_CALC_DATA_FILE, true, "Optional: file containing computed data per sample");
        options.addOption(BA_SPECIFIC_CANCER, true, "Optional: Only process this cancer type");
    }

    public boolean initialise(GenericDataCollection collection, final CommandLine cmd)
    {
        mDataCollection = collection;
        mOutputFileId = cmd.getOptionValue(OUTPUT_FILE_ID);
        mOutputDir = cmd.getOptionValue(OUTPUT_DIR);
        mSpecificCancer = cmd.getOptionValue(BA_SPECIFIC_CANCER, "");

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

            sample.setSampleName(mDataCollection.getFieldNames().get(sampleId));

            final List<String> extSampleData = getSampleExtData(sample.getSampleName());
            if(extSampleData != null)
            {
                sample.setCancerType(extSampleData.get(COL_CANCER_TYPE));
                sample.setCategoryData(extSampleData);
            }

            if(!mSpecificCancer.isEmpty() && !sample.getCancerType().equals(mSpecificCancer))
            {
                sample.setExcluded(true);
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
        calcCountsNoise();
        collectElevatedSampleBuckets();
        perfCounter.stop();

        // back-ground groups using expected counts
        perfCounter.start("BackgroundGroups");
        formBackgroundBucketGroups();
        // logBucketGroups(false, true);
        perfCounter.stop();

        int maxRuns = 25;

        for(int runId = 0; runId < maxRuns; ++runId)
        {
            perfCounter.start(String.format("FindBucketGroup run %d", runId));

            mBucketGroups.clear();

            LOGGER.debug("starting run({}) to find next top bucket group", runId);

            // search for precise bucket groups
            formExactBucketGroups();

            // and then the best sub-groups
            formBucketGroupsFromSamplePairs();

            double prevAllocCount = mAllocatedCount;

            populateTopBucketGroups();

            // allocated elevated counts to the best group
            BucketGroup nextBestGroup = allocateTopBucketGroup();

            if(nextBestGroup != null)
            {
                logOverallStats();
                nextBestGroup.setSelected(true);
                mFinalBucketGroups.add(nextBestGroup);
            }

            perfCounter.stop();

            mBucketGroups.clear();

            // temp for logging only
            mBucketGroups.addAll(mFinalBucketGroups);
            analyseGroupsVsExtData(true, false);
            logBucketGroups(true);

            if (nextBestGroup == null)
                break;

            if((mAllocatedCount - prevAllocCount) / mElevatedCount < 0.001)
                break;
        }

        perfCounter.start("FinalFit");
        fitAllSamples();
        logBucketGroups(true);
        perfCounter.stop();

        // start by looking for near-exact matches in buckets across samples, forming bucket groups
        // formExactBucketGroups();
        analyseGroupsVsExtData(true, true); // just to see the splits
        //logBucketGroups(true, true);

        perfCounter.start("AnalyseResults");
        // analyseGroupsVsExtData(true, true);
        logSampleResults();
        logWorstAllocatedSamples();
        perfCounter.stop();

//        createSignatures();
//        compareSignatures();

        writeBucketGroups();

        writeSampleData();
        writeBackgroundSigs();
        writeSampleMatrixData(mElevatedCounts, "_ba_elevated_counts.csv");
        writeSampleCalcData();

        perfCounter.logStats();

        finalise();
    }

    private void finalise()
    {
        try
        {
            if(mBucketGroupFileWriter != null)
                mBucketGroupFileWriter.close();
        }
        catch(IOException e)
        {
            LOGGER.error("failed to close bucket groups file");
        }
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

        // double bucketMedRange = 0.005;

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

        double[][] brData = mBucketMedianRatios.getData();
        double[][] probData = mBucketProbs.getData();
        double[][] bgData = mBackgroundCounts.getData();
        double[][] elevData = mElevatedCounts.getData();
        final double[][] scData = mSampleCounts.getData();

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

                    int backgroundCount = (int)round(bucketMedianRatio * bgAllocation);

                    // expData[i][j] = expectedCount;
                    // if(expectedCount < backgroundCount)
                    //    backgroundCount = expectedCount;

                    if(sbCount < backgroundCount) // must be capped by the actual count
                        backgroundCount = sbCount;

                    bgData[i][j] = backgroundCount;

                    if (backgroundCount > 0)
                    {
                        brData[i][j] = sbCount / (double)backgroundCount;
                    }

                    int elevatedCount = sbCount - backgroundCount;

                    // compute a range for Poisson noise around this elevated count
                    if (elevatedCount > 0)
                    {
                        elevData[i][j] = elevatedCount;
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

        // if configured, switch all background counts to elevated
        if(mNoBackgroundCounts)
        {
            for (int i = 0; i < mBucketCount; ++i)
            {
                for (int j = 0; j < mSampleCount; ++j)
                {
                    elevData[i][j] += bgData[i][j];
                    bgData[i][j] = 0;
                }
            }
        }

        mBackgroundCounts.cacheTranspose();
        mElevatedCounts.cacheTranspose();

        mBackgroundCount = mBackgroundCounts.sum();
        mElevatedCount = mElevatedCounts.sum();

        LOGGER.debug(String.format("total counts: background(%s perc=%.3f) elevated(%s perc=%.3f) of total(%s)",
                sizeToStr(mBackgroundCount), mBackgroundCount/mTotalCount,
                sizeToStr(mElevatedCount), mElevatedCount/mTotalCount, sizeToStr(mTotalCount)));

        for (int i = 0; i < probFrequencies.length-1; ++i)
        {
            LOGGER.debug(String.format("probability(1e-%d) freq(%d) percOfTotal(%.4f)",
                    i, probFrequencies[i], probFrequencies[i]/(double)gridSize));
        }

        LOGGER.debug(String.format("probability(zero) freq(%d) percOfTotal(%.4f)",
                probFrequencies[zeroProbIndex], probFrequencies[zeroProbIndex]/(double)gridSize));
    }

    private void calcCountsNoise()
    {
        mPermittedElevRange = new NmfMatrix(mBucketCount, mSampleCount);
        mPermittedBgRange = new NmfMatrix(mBucketCount, mSampleCount);

        if(!mApplyNoise)
            return;

        double[][] bgData = mBackgroundCounts.getData();
        double[][] elevData = mElevatedCounts.getData();
        double[][] permBgRangeData = mPermittedBgRange.getData();
        double[][] permElevRangeData = mPermittedElevRange.getData();

        Map<Integer, Integer> rangeMap = new HashMap(); // cache range values

        for (int i = 0; i < mBucketCount; ++i)
        {
            for (int j = 0; j < mSampleCount; ++j)
            {
                int backgroundCount = (int)bgData[i][j];

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
                }

                int elevatedCount = (int)elevData[i][j];

                // compute a range for Poisson noise around this elevated count
                if (elevatedCount > 0)
                {
                    elevData[i][j] = elevatedCount;

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
            }
        }

        mPermittedElevRange.cacheTranspose();
        mPermittedBgRange.cacheTranspose();
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
            sample.setElevatedBucketCounts(mElevatedCounts.getCol(i), mPermittedElevRange.getCol(i));

            if (!bucketList.isEmpty())
            {
                sample.setElevatedBuckets(bucketList);
                mAllSampleBucketGroups.put(i, bucketList);

//                if(sample.getElevatedCount() == 0)
//                {
//                    LOGGER.debug(String.format("sample(%d) has %d elevated bucket but zero elevated count, total(%s bg=%s)",
//                            i, bucketList.size(), sizeToStr(sample.getTotalCount()), sizeToStr(sumVector(mBackgroundCounts.getCol(i)))));
//                }
            }
        }

        // recalc totals just for the applicable cancer type
        if(!mSpecificCancer.isEmpty())
        {
            mTotalCount = 0;
            mBackgroundCount = 0;
            mElevatedCount = 0;

            // calculate manually
            for(final SampleData sample : mSampleData)
            {
                if(sample.isExcluded())
                    continue;

                mTotalCount += sample.getTotalCount();
                mBackgroundCount += sample.getTotalCount() - sample.getElevatedCount();
                mElevatedCount += sample.getElevatedCount();
            }
        }

        LOGGER.debug(String.format("samples with elevated buckets: count(%d perc=%.2f), buckets(%d perc=%.3f)",
                mAllSampleBucketGroups.size(), mAllSampleBucketGroups.size()/(double)mSampleCount,
                totalCount, totalCount/(double)(mBucketCount*mSampleCount)));
    }

    private void formExactBucketGroups()
    {
        // forms groups out of samples with similarly elevated buckets
        int groupsCreated = 0;

        for (int samIndex1 = 0; samIndex1 < mSampleCount; ++samIndex1)
        {
            SampleData sample1 = mSampleData.get(samIndex1);

            if(sample1.isExcluded())
                continue;

            if(!mReassessSamples.isEmpty() && !mReassessSamples.contains(samIndex1))
                continue;

            final List<Integer> bl1 = sample1.getUnallocBuckets();

            if (bl1.isEmpty())
                continue;

            if(sample1.getUnallocPercent() < minAllocPercent(sample1))
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

                    double[] sc1 = extractBucketCountSubset(sample1, combinedBuckets);
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
                    bucketGroup.addSample(samIndex1, sc1);

                    if(mLogVerbose)
                    {
                        LOGGER.debug(String.format("bg(%d) added sample(%d) with buckets(grp=%d sam=%d match=%d) css(%.4f prob=%s) totalSamples(%d)",
                                bucketGroup.getId(), samIndex1, bucketGroup.getBucketIds().size(), bl1.size(),
                                commonBuckets.size(), bcCss, withinProb ? "pass" : "fail", bucketGroup.getSampleIds().size()));
                    }

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

                if(sample2.isExcluded())
                    continue;;

                final List<Integer> bl2 = sample2.getUnallocBuckets();

                if (bl2.isEmpty())
                    continue;

                if(sample2.getUnallocPercent() < minAllocPercent(sample2))
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

                // work out the proportion of elevated counts are covered by these shared elevated buckets
                double sharedElevatedCount = 0;
                final double[] sam1ElevCounts = sample1.getUnallocBucketCounts();
                final double[] sam2ElevCounts = sample2.getUnallocBucketCounts();
                for(Integer bucket : commonBuckets)
                {
                    sharedElevatedCount += sam1ElevCounts[bucket] + sam2ElevCounts[bucket];
                }

                double sharedElevCountPerc = sharedElevatedCount / (sample1.getUnallocatedCount() + sample2.getUnallocatedCount());
                boolean closeBucketMatch = (minMatch >= mHighRequiredMatch) || sharedElevCountPerc >= mHighRequiredMatch;

                if (closeBucketMatch)
                {
                    double[] sc1 = extractBucketCountSubset(sample1, commonBuckets);
                    double[] sc2 = extractBucketCountSubset(sample2, commonBuckets);
                    double bcCss = calcSharedCSS(sc1, sc2);

                    if(bcCss >= mHighCssThreshold)
                    {
                        BucketGroup bucketGroup = new BucketGroup(mNextBucketId++);
                        bucketGroup.setTag("exact");
                        bucketGroup.addBuckets(commonBuckets);

                        bucketGroup.addSample(samIndex1, sc1);
                        bucketGroup.addSample(samIndex2, sc2);

                        if(mLogVerbose)
                        {
                            LOGGER.debug(String.format("added bg(%d) samples(%d and %d) with buckets(s1=%d s2=%d match=%d) sharedElevPerc(%.3f) css(%.4f)",
                                    bucketGroup.getId(), samIndex1, samIndex2, bl1.size(), bl2.size(), commonBucketCount, sharedElevCountPerc, bcCss));
                        }

                        refineBucketGroupBuckets(bucketGroup);

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
            LOGGER.debug("no exact-match bucket groups created");
            return;
        }

        LOGGER.debug("{} exact-match bucket groups created", groupsCreated);

//        LOGGER.debug("bucket groups created({}) additions({}) total({}), unallocated elevated samples({})",
//                groupsCreated, groupsAdjusted, mBucketGroups.size(), mSampleCount - (groupsCreated * 2) - groupsAdjusted);
    }

    private void formBucketGroupsFromSamplePairs()
    {
        // create unique groups from any subset of 2 samples' buckets if they are a close match
        // samples aren't necessarily allocated, only groups are made for the time-being
        LOGGER.debug("forming bucket groups from sample-pair reduced buckets");

        int groupsCreated = 0;

        for (int samIndex1 = 0; samIndex1 < mSampleCount; ++samIndex1)
        {
            SampleData sample1 = mSampleData.get(samIndex1);

            if(sample1.isExcluded())
                continue;

            if(!mReassessSamples.isEmpty() && !mReassessSamples.contains(samIndex1))
                continue;

            if(mSampleWatchList.contains(samIndex1))
            {
                //LOGGER.debug("spec sample");
            }

            double reqSam1AllocPercent = minAllocPercent(sample1);

            if(sample1.getUnallocPercent() < reqSam1AllocPercent)
                continue;

            final List<Integer> bl1 = sample1.getUnallocBuckets();

            if (bl1.isEmpty())
                continue;

            // record the top matching other sample - this will be used to create the top-allocating group
            double maxAllocaTotal = 0;
            int maxOtherSample = -1;
            List<Double> maxCombinedCounts = Lists.newArrayList();
            List<Integer> maxSharedBuckets = Lists.newArrayList();

            for (int samIndex2 = samIndex1 + 1; samIndex2 < mSampleCount; ++samIndex2)
            {
                SampleData sample2 = mSampleData.get(samIndex2);

                if(sample2.isExcluded())
                    continue;;

                double reqSam2AllocPercent = minAllocPercent(sample2);

                if(sample2.getUnallocPercent() < reqSam2AllocPercent)
                    continue;

                final List<Integer> bl2 = sample2.getUnallocBuckets();

                if (bl2.isEmpty())
                    continue;

                List<Integer> commonBuckets = getMatchingList(bl1, bl2);
                int commonBucketCount = commonBuckets.size();

                if (commonBucketCount < 2)
                    continue;

                double[] sc1 = extractBucketCountSubset(sample1, commonBuckets);
                double[] sc2 = extractBucketCountSubset(sample2, commonBuckets);

                final double[] sam1ElevCounts = sample1.getUnallocBucketCounts();
                final double[] sam2ElevCounts = sample2.getUnallocBucketCounts();

                double allBucketsTotal = sample1.getElevatedCount() + sample2.getElevatedCount();
                double sam1ElevTotal = sumVector(sc1);
                double sam2ElevTotal = sumVector(sc2);
                double elevatedTotal = sam1ElevTotal + sam2ElevTotal;

                if(sam1ElevTotal/sample1.getElevatedCount() < reqSam1AllocPercent || sam2ElevTotal/sample2.getElevatedCount() < reqSam2AllocPercent)
                    continue;

                // only create groups form overlapping buckets which contribute at least half the samples' mutation load
//                if(elevatedTotal < allBucketsTotal * 0.5)
//                    continue;

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
                    double[] cssResults = new double[commonBuckets.size()];

                    for (int i = 0; i < commonBuckets.size(); ++i)
                    {
                        Integer testBucket = commonBuckets.get(i);

                        sc1[testBucket] = 0;
                        sc2[testBucket] = 0;

                        // run CSS on this reduced set of buckets
                        cssResults[i] = calcSharedCSS(sc1, sc2);

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
                        bcCss = calcSharedCSS(sc1, sc2);

                        removedBuckets.add(testBucket);

                        if (bcCss >= mHighCssThreshold)
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

                double newAllocTotal = elevatedTotal * pow(bcCss, 2);

                if(newAllocTotal <= maxAllocaTotal)
                    continue;

                maxAllocaTotal = newAllocTotal;
                maxOtherSample = samIndex2;
                maxCombinedCounts.clear();;
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
            BucketGroup bucketGroup = new BucketGroup(mNextBucketId++);
            bucketGroup.setTag("subset");

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

            if(mLogVerbose)
            {
                LOGGER.debug(String.format("added bg(%d) samples(%d and %d) with buckets(%d) allocCalcTotal(%s)",
                        bucketGroup.getId(), samIndex1, maxOtherSample, maxSharedBuckets.size(), sizeToStr(maxAllocaTotal)));
            }

            // refineBucketGroupBuckets(bucketGroup);

            mBucketGroups.add(bucketGroup);
            ++groupsCreated;
        }

        if(mBucketGroups.isEmpty())
        {
            LOGGER.debug("no sample-pair subset bucket groups created");
            return;
        }

        LOGGER.debug("sample-pair subset created({}))", groupsCreated);
    }

    private void refineBucketGroupBuckets(BucketGroup bucketGroup)
    {
        if(bucketGroup.getBucketIds().size() <= 5)
            return;

        double minPercentOfMax = 0.001;
        double minRatioTotal = 0.95;
        double bucketRatioTotal = 0;
        double maxBucketRatio = 0;

        double[] bucketRatios = bucketGroup.getBucketRatios();

        List<Integer> sortedRatioIndices = getSortedVectorIndices(bucketRatios, false);
        List<Integer> newBucketsList = Lists.newArrayList();

        for(Integer bucketId : sortedRatioIndices)
        {
            double bucketRatio = bucketRatios[bucketId];

            if(maxBucketRatio == 0)
                maxBucketRatio = bucketRatio;

            if(bucketRatio == 0)
                break;

            if(bucketRatioTotal >= minRatioTotal && bucketRatio < minPercentOfMax * maxBucketRatio)
                break;

            bucketRatioTotal += bucketRatio;
            newBucketsList.add(bucketId);
        }

        if(newBucketsList.size() == bucketGroup.getBucketIds().size())
            return;

//        LOGGER.debug(String.format("bg(%d) refined buckets(%d -> %d) ratios(total=%.3f max=%.3f avg=%.3f)",
//                bucketGroup.getId(), bucketGroup.getBucketIds().size(), newBucketsList.size(),
//                bucketRatioTotal, maxBucketRatio, bucketRatioTotal/newBucketsList.size()));

        bucketGroup.reduceToBucketSet(newBucketsList);
    }

    private void populateTopBucketGroups()
    {
        mTopAllocBucketGroups.clear();

        LOGGER.debug("finding top potential bucket group from count({})", mBucketGroups.size());

        SigContributionOptimiser sigOptim = new SigContributionOptimiser(mBucketCount, false, MIN_GROUP_ALLOC_PERCENT_LOWER, SAMPLE_ALLOCATED_PERCENT, true);

        // first clear all existing allocations of samples to groups and vice versa
        for (BucketGroup bucketGroup : mBucketGroups)
        {
            bucketGroup.clearSamples();
            bucketGroup.resetPotentialAllocation();

            double[] bgRatios = new double[mBucketCount];
            copyVector(bucketGroup.getBucketRatios(), bgRatios);

            final List<Integer> groupBuckets = bucketGroup.getBucketIds();

            for (int sampleId = 0; sampleId < mSampleCount; ++sampleId)
            {
                final SampleData sample = mSampleData.get(sampleId);

                if(sample.isExcluded())
                    continue;

                if(mSampleWatchList.contains(sampleId))
                {
                    //LOGGER.debug("spec sample");
                }

                final List<Integer> samBuckets = sample.getElevatedBuckets();

                double reqAllocPercent = minAllocPercent(sample);

                if (samBuckets.isEmpty() || sample.getUnallocPercent() < reqAllocPercent)
                    continue;

                final List<Integer> commonBuckets = getMatchingList(groupBuckets, samBuckets);

                // ensure sample's buckets are covered by the group, even if the sample has more
                double groupPerc = commonBuckets.size() / (double) groupBuckets.size();
                if (groupPerc < mHighRequiredMatch)
                {
                    // LOGGER.debug(String.format("sample(%d) excluded on low groupPerc(%.3f)", sampleId, groupPerc));
                    continue;
                }

                double[] allocCounts = sample.getPotentialUnallocCounts(bgRatios, groupBuckets);

                int allocCountTotal = (int)round(sumVector(allocCounts));
                double allocPercent = allocCountTotal / sample.getElevatedCount();

                if (allocPercent < reqAllocPercent)
                {
                    // see if a fit with sig along with all the other allocated one for this sample would then meet the min % threshold
                    if (!sample.getElevBucketGroups().isEmpty())
                    {
                        List<double[]> ratiosCollection = Lists.newArrayList();

                        for (final BucketGroup samGroup : sample.getElevBucketGroups())
                        {
                            ratiosCollection.add(samGroup.getBucketRatios());
                        }

                        ratiosCollection.add(bgRatios);

                        double[] prevContribs = new double[ratiosCollection.size()];
                        int candidateSigIndex = prevContribs.length - 1;

                        sigOptim.initialise(sample.Id, sample.getElevatedBucketCounts(), sample.getCountRanges(), ratiosCollection, prevContribs);
                        sigOptim.setTargetSig(candidateSigIndex);
                        boolean validCalc = sigOptim.fitToSample(SAMPLE_ALLOCATED_PERCENT, reqAllocPercent);

                        if (!validCalc) // couldn't reach the required percent for this canidate sig
                            continue;

                        allocPercent = sigOptim.getContribs()[candidateSigIndex] / sample.getElevatedCount();
                    }

                    if (allocPercent < reqAllocPercent)
                        continue;
                }

                bucketGroup.addPotentialAllocation(allocCountTotal);
                bucketGroup.addPotentialAdjAllocation(allocCountTotal * allocPercent);
                bucketGroup.addSample(sampleId, allocCounts);
            }

            if(bucketGroup.getPotentialAllocation() == 0)
                continue;

            // store in order since the top one(s) will be allocated first
            int bgIndex = 0;
            while(bgIndex < mTopAllocBucketGroups.size())
            {
                if(bucketGroup.getPotentialAdjAllocation() >= mTopAllocBucketGroups.get(bgIndex).getPotentialAdjAllocation())
                    break;

                ++bgIndex;
            }

            if(bgIndex < 50)
                mTopAllocBucketGroups.add(bgIndex, bucketGroup);
        }

        if(mTopAllocBucketGroups.isEmpty())
        {
            LOGGER.warn("no max potential group found");
        }

        LOGGER.debug("found {} top bucket groups", mTopAllocBucketGroups.size());
    }

    private BucketGroup allocateTopBucketGroup()
    {
        if(mTopAllocBucketGroups.isEmpty())
            return null;

        mReassessSamples.clear();

        BucketGroup topBucketGroup = mTopAllocBucketGroups.get(0);

        LOGGER.debug(String.format("top bg(%d) with buckets(%d) samples(%d) potential allocation(%s adj=%s)",
                topBucketGroup.getId(), topBucketGroup.getBucketIds().size(), topBucketGroup.getSampleIds().size(),
                sizeToStr(topBucketGroup.getPotentialAllocation()), sizeToStr(topBucketGroup.getPotentialAdjAllocation())));

        final double[] topBucketRatios = topBucketGroup.getBucketRatios();
        double[] avgBucketRatios = new double[mBucketCount];
        copyVector(topBucketRatios, avgBucketRatios);

        List<Integer> topSampleIds = Lists.newArrayList();
        topSampleIds.addAll(topBucketGroup.getSampleIds());

        List<Double> sampleAllocTotals = Lists.newArrayList();
        sampleAllocTotals.addAll(topBucketGroup.getSampleCountTotals());

        // walk down list, recalculating the ratios, expanding their range and comparing to groups further down
        int maxComparisons = 50;
        double topPotentialAlloc = topBucketGroup.getPotentialAllocation();
        List<Integer> candiateAdditionalBuckets = Lists.newArrayList();

        for(int bgIndex = 1; bgIndex < min(mTopAllocBucketGroups.size(), maxComparisons); ++bgIndex)
        {
            BucketGroup bucketGroup = mTopAllocBucketGroups.get(bgIndex);

            // recompute ratios and widen
            double[] bucketRatios = bucketGroup.getBucketRatios();

            double groupCss = calcCSS(topBucketRatios, bucketRatios);
            if(groupCss < mHighCssThreshold)
                continue;

            double groupToTopRatio = bucketGroup.getPotentialAllocation() / topPotentialAlloc;
            vectorMultiply(bucketRatios, groupToTopRatio);

            for(Integer bucket : topBucketGroup.getBucketIds())
            {
                avgBucketRatios[bucket] += bucketRatios[bucket];
            }

            List<Integer> missingBuckets = getDiffList(bucketGroup.getBucketIds(), topBucketGroup.getBucketIds());

            for (Integer bucket : missingBuckets)
            {
                if(!candiateAdditionalBuckets.contains(bucket))
                    candiateAdditionalBuckets.add(bucket);
            }

            // merge in any difference in samples
            final List<Integer> bgSamples = bucketGroup.getSampleIds();
            final List<Double> bgSampleTotals = bucketGroup.getSampleCountTotals();

            int samplesAdded = 0;
            for(int samIndex = 0; samIndex < bgSamples.size(); ++samIndex)
            {
                Integer sampleId = bgSamples.get(samIndex);

                if (!topSampleIds.contains(sampleId))
                {
                    topSampleIds.add(sampleId);
                    sampleAllocTotals.add(bgSampleTotals.get(samIndex));
                    ++samplesAdded;
                }
            }

            LOGGER.debug(String.format("top bg(%d) merged with bg(%d) alloc(%s adj=%s) css(%.4f) samples(bg1=%d bg2=%d added=%d)",
                    topBucketGroup.getId(), bucketGroup.getId(),
                    sizeToStr(bucketGroup.getPotentialAllocation()), sizeToStr(bucketGroup.getPotentialAdjAllocation()),
                    groupCss, topSampleIds.size(), bgSamples.size(), samplesAdded));

            // remove this from the candidate set so it doesn't impact the skipped-sample logic
            for(final BucketGroup candidateBg : mBucketGroups)
            {
                if(candidateBg == bucketGroup)
                {
                    mBucketGroups.remove(candidateBg);
                    break;
                }
            }
        }

        /*
        // test out any candidate additional buckets
        if(!candiateAdditionalBuckets.isEmpty())
        {
            double[] sampleAllocTotalsArray = listToArray(sampleAllocTotals);

            double[] newBucketRatios = new double[mBucketCount];
            copyVector(topBucketRatios, newBucketRatios);
            List<Integer> newBuckets = Lists.newArrayList();
            newBuckets.addAll(topBucketGroup.getBucketIds());


            for (Integer bucket : candiateAdditionalBuckets)
            {
                double[] sbCounts = new double[topSampleIds.size()];

                for (int samIndex = 0; samIndex < topSampleIds.size(); ++samIndex)
                {
                    final SampleData sample = mSampleData.get(topSampleIds.get(samIndex));
                    sbCounts[samIndex] = sample.getUnallocBucketCounts()[bucket];
                }

                double bucketCss = calcCSS(sbCounts, sampleAllocTotalsArray);

//                if(bucketCss <= mHighCssThreshold)
//                    continue;

                // convert to a ratio and then test against each bucket's unallocated counts
                // since we don't want it to be a limited factor
                double lsRatio = calcLinearLeastSquares(sampleAllocTotalsArray, sbCounts);
                double minRatio = calcMinPositiveRatio(sampleAllocTotalsArray, sbCounts);
                double rawRatio = sumVector(sbCounts) / sumVector(sampleAllocTotalsArray);

                newBucketRatios[bucket] = lsRatio;
                convertToPercentages(newBucketRatios);
                newBuckets.add(bucket);

                double totalAllocChange = 0;

                // now test against every sample to see if adding this new bucket would lower the potential allocation
                for (int samIndex = 0; samIndex < topSampleIds.size(); ++samIndex)
                {
                    final SampleData sample = mSampleData.get(topSampleIds.get(samIndex));
                    double prevAlloc = sampleAllocTotalsArray[samIndex];

                    double[] newaAllocCounts = sample.getPotentialAllocation(newBucketRatios, newBuckets);
                    double newAllocTotal = sumVector(newaAllocCounts);

                    totalAllocChange += (newAllocTotal - prevAlloc);
                }

                newBucketRatios[bucket] = 0;
                newBuckets.remove(newBuckets.size() - 1);

                LOGGER.debug(String.format("candidate bucket(%d) css(%.4f) ratio(raw=%.3f leastSq=%.3f min=%.3f) netAllocChange(%s)",
                        bucket, bucketCss, rawRatio, lsRatio, minRatio, sizeToStr(totalAllocChange)));

                if(totalAllocChange > 0)
                {
                    LOGGER.debug(String.format("top bg(%d) adding bucket(%d) with net positive allocation(%s)",
                            topBucketGroup.getId(), bucket, sizeToStr(totalAllocChange)));

                    //topBucketGroup.addBucket(bucket, false);
                    // avgBucketRatios[bucket] = lsRatio;
                }
            }
        }
        */

        topBucketGroup.clearSamples();

        // convert back to percentages
        convertToPercentages(avgBucketRatios);
        topBucketGroup.setBucketRatios(avgBucketRatios);

        calcBucketRatioRanges(topBucketGroup, topSampleIds, sampleAllocTotals);

        final List<Integer> groupBuckets = topBucketGroup.getBucketIds();
        final List<Double> bucketRatioRanges = topBucketGroup.getBucketRatioRanges();

        // now allocate samples to this top group
        double totalAlloc = 0;

        List<Integer> skippedSamples = Lists.newArrayList();

        for(int samIndex = 0; samIndex < topSampleIds.size(); ++samIndex)
        {
            Integer sampleId = topSampleIds.get(samIndex);
            double proposedAlloc = sampleAllocTotals.get(samIndex);

            final SampleData sample = mSampleData.get(sampleId);

            if(mSampleWatchList.contains(sampleId))
            {
                //LOGGER.debug("spec sample");
            }

            double[] sampleCountAllocations = sample.getPotentialAllocation(avgBucketRatios, groupBuckets, bucketRatioRanges);

            double newAllocTotal = sumVector(sampleCountAllocations);

            double maxOtherGroupAlloc = getSampleMaxAllocation(sample, topBucketGroup, true);

            if(maxOtherGroupAlloc > SKIP_ALLOC_FACTOR * newAllocTotal)
            {
                LOGGER.debug(String.format("sample(%d) skipped bg(%d) with better allocation(this=%s other=%s)",
                        sampleId, topBucketGroup.getId(), sizeToStr(newAllocTotal), sizeToStr(maxOtherGroupAlloc)));

                skippedSamples.add(sampleId);
                continue;
            }

            double reqAllocPercent = minAllocPercent(sample);

            // apply to the remaining unallocated elevated counts for this sample
            double prevAllocPerc = sample.getAllocPercent();
            double actualAlloc = sample.allocateBucketCounts(sampleCountAllocations, reqAllocPercent);
            double allocPerc = actualAlloc / sample.getElevatedCount();

            if(allocPerc >= reqAllocPercent)
            {
                topBucketGroup.addSample(sampleId, sampleCountAllocations, false);
                sample.addElevBucketGroup(topBucketGroup, allocPerc);
                totalAlloc += actualAlloc;

                LOGGER.debug(String.format("sample(%d) added to bg(%d) buckets(grp=%d sam=%d unalloc=%d) count(orig=%s rcRng=%s act=%s of %s) allocatedPerc(%.3f -> %.3f) groupCount(%d)",
                        sampleId, topBucketGroup.getId(), groupBuckets.size(), sample.getElevatedBuckets().size(), sample.getUnallocBuckets().size(),
                        sizeToStr(proposedAlloc), sizeToStr(newAllocTotal), sizeToStr(actualAlloc), sizeToStr(sample.getElevatedCount()),
                        prevAllocPerc, sample.getAllocPercent(), sample.getElevBucketGroups().size()));
            }
        }

        LOGGER.debug(String.format("new top bg(%d) added %d samples, totalAllocatedCount(%s)",
                topBucketGroup.getId(), topBucketGroup.getSampleIds().size(), sizeToStr(totalAlloc)));

        // recheck any sample previously skipped against the fixed groups so see if they can now be allocated
        checkSkippedSamples(skippedSamples, topBucketGroup.getSampleIds());

        // for now clear all top groups to force a reassessment, since after allocation of elevated counts
        // to the top bucket, the similarities could increase across more buckets and samples
        // purgeSimilarBucketGroups(topBucketGroup);
        mTopAllocBucketGroups.clear();

        if(topBucketGroup.getSampleIds().isEmpty())
            return null;

        mReassessSamples.addAll(topBucketGroup.getSampleIds());

        return topBucketGroup;
    }

    private double getSampleMaxAllocation(final SampleData sample, final BucketGroup excludeGroup, boolean checkFinalGroups)
    {
        double maxAlloc = 0;

        // test against the candidate list of bucket groups
        BucketGroup bestGroup = null;

        for(final BucketGroup bucketGroup : mBucketGroups)
        {
            if(bucketGroup == excludeGroup)
                continue;

            for (int samIndex = 0; samIndex < bucketGroup.getSampleIds().size(); ++samIndex)
            {
                int sampleId = bucketGroup.getSampleIds().get(samIndex);

                if (sampleId == sample.Id)
                {
                    double alloc = bucketGroup.getSampleCountTotals().get(samIndex);
                    if (alloc > maxAlloc)
                    {
                        maxAlloc = alloc;
                        bestGroup = bucketGroup;
                    }
                }
            }
        }

        if(checkFinalGroups)
        {
            // also test against the final set of bucket groups
            for (final BucketGroup bucketGroup : mFinalBucketGroups)
            {
                if(bucketGroup.hasSample(sample.Id))
                    continue;

                double[] allocCounts = sample.getPotentialUnallocCounts(bucketGroup.getBucketRatios(), bucketGroup.getBucketIds());

                double alloc = sumVector(allocCounts);

                if (alloc > maxAlloc)
                {
                    maxAlloc = alloc;
                }
            }
        }

//        if(bestGroup != null)
//        {
//            LOGGER.debug(String.format("sample(%d) best bg(%d) with maxAlloc(%s)",
//                    sample.Id, bestGroup.getId(), sizeToStr(maxAlloc)));
//        }

        return maxAlloc;
    }

    private void checkSkippedSamples(List<Integer> newSkippedSamples, List<Integer> newAddedSamples)
    {
        // for samples previously skipped, check if can now be allocated to one or more of the final bucket groups
        if(mFinalBucketGroups.isEmpty())
            return;

        if(newSkippedSamples.isEmpty() && mSkippedSamples.isEmpty())
            return;

        int initSkippedSamples = mSkippedSamples.size();

        // first remove the ones just added to a group
        for(Integer sampleId : newAddedSamples)
        {
            if(mSkippedSamples.contains(sampleId))
                mSkippedSamples.remove(sampleId);
        }

        int samIndex = 0;
        while(samIndex < mSkippedSamples.size())
        {
            Integer sampleId = mSkippedSamples.get(samIndex);

            if(newSkippedSamples.contains(sampleId))
            {
                ++samIndex;
                continue;
            }

            SampleData sample = mSampleData.get(sampleId);

            double maxOtherGroupAlloc = getSampleMaxAllocation(sample, null, false);
            double maxOtherGroupAllocPerc = maxOtherGroupAlloc / sample.getElevatedCount();

            double reqAllocPercent = minAllocPercent(sample);

            boolean wasAllocated = false;

            for(final BucketGroup bucketGroup : mFinalBucketGroups)
            {
                if(bucketGroup.hasSample(sampleId))
                    continue;

                double[] allocCounts = sample.getPotentialAllocation(bucketGroup.getBucketRatios(), bucketGroup.getBucketIds(), bucketGroup.getBucketRatioRanges());

                double proposedAllocTotal = sumVector(allocCounts);
                double proposedAllocPerc = proposedAllocTotal / sample.getElevatedCount();

                if(proposedAllocPerc < reqAllocPercent)
                    continue;

                // continue to hold out if there is a better candidate group yet to come
                if(maxOtherGroupAlloc > SKIP_ALLOC_FACTOR * proposedAllocTotal)
                    continue;

                // apply to the remaining unallocated elevated counts for this sample
                double prevAllocPerc = sample.getAllocPercent();
                double actualAlloc = sample.allocateBucketCounts(allocCounts, reqAllocPercent);
                double allocPerc = actualAlloc / sample.getElevatedCount();

                if(allocPerc >= reqAllocPercent)
                {
                    bucketGroup.addSample(sampleId, allocCounts, false);
                    sample.addElevBucketGroup(bucketGroup, allocPerc);
                    wasAllocated = true;
                    mReassessSamples.add(sampleId);

                    LOGGER.debug(String.format("sample(%d) skipped now added to bg(%d) buckets(grp=%d sam=%d unalloc=%d) count(prop=%s act=%s of %s) allocatedPerc(%.3f -> %.3f) groupCount(%d)",
                            sampleId, bucketGroup.getId(), bucketGroup.getBucketIds().size(), sample.getElevatedBuckets().size(), sample.getUnallocBuckets().size(),
                            sizeToStr(proposedAllocTotal), sizeToStr(actualAlloc), sizeToStr(sample.getElevatedCount()),
                            prevAllocPerc, sample.getAllocPercent(), sample.getElevBucketGroups().size()));
                }
            }

            if(!wasAllocated)
            {
                if(maxOtherGroupAllocPerc >= reqAllocPercent)
                {
//                    LOGGER.debug(String.format("sample(%d) skipped again with better allocation(%s perc=%.3f)",
//                            sampleId, sizeToStr(maxOtherGroupAlloc), maxOtherGroupAllocPerc));

                    ++samIndex;
                    continue;
                }
                else
                {
                    LOGGER.debug("sample({}) skipped previously but not allocated later", sampleId);
                }
            }

            mSkippedSamples.remove(samIndex);
        }

        // finally add in the newly skipped samples
        for(Integer sampleId : newSkippedSamples)
        {
            if(!mSkippedSamples.contains(sampleId))
                mSkippedSamples.add(sampleId);
        }

        if(initSkippedSamples > mSkippedSamples.size())
        {
            LOGGER.debug("skipped samples added({}), overall({} -> {})", newSkippedSamples.size(), initSkippedSamples, mSkippedSamples.size());
        }
    }

    private void fitAllSamples()
    {
        LOGGER.debug("applying {} final bucket groups to all samples from scratch", mFinalBucketGroups.size());

        for(BucketGroup bucketGroup : mFinalBucketGroups)
        {
            bucketGroup.clearSamples();
        }

        double reqAllocPercent = MIN_GROUP_ALLOC_PERCENT_LOWER;

        SigContributionOptimiser sigOptim = new SigContributionOptimiser(mBucketCount, true, MIN_GROUP_ALLOC_PERCENT_LOWER, SAMPLE_ALLOCATED_PERCENT, true);

        if(mNoBackgroundCounts)
        {
            LOGGER.debug("including {} background group(s)", mBackgroundGroups.size());

            List<BucketGroup> elevatedGroups = Lists.newArrayList();
            elevatedGroups.addAll(mFinalBucketGroups);

            // put the BG groups first
            mFinalBucketGroups.clear();
            mFinalBucketGroups.addAll(mBackgroundGroups);
            mFinalBucketGroups.addAll(elevatedGroups);
        }

        for (SampleData sample : mSampleData)
        {
            if (sample.isExcluded() || sample.getElevatedCount() == 0)
                continue;

            double prevTotalAllocPerc = sample.getAllocPercent();
            int prevGroupCount = sample.getElevBucketGroups().size();

            sample.clearAllocations();

            double sampleCount = sample.getElevatedCount();

            List<Integer> bgIndexList = Lists.newArrayList();
            List<Double> potentialAllocTotals = Lists.newArrayList();
            List<double[]> potentialAllocCounts = Lists.newArrayList();

            if(mSampleWatchList.contains(sample.Id))
            {
                 LOGGER.debug("spec sample");
            }

            for(int bgIndex = 0; bgIndex < mFinalBucketGroups.size(); ++bgIndex)
            {
                BucketGroup bucketGroup = mFinalBucketGroups.get(bgIndex);

                if(mBackgroundGroups.contains(bucketGroup))
                {
                    // assign if applicable by cancer type
                    if(sample.getBackgroundGroup() != bucketGroup)
                        continue;
                }

                // re-test with all elevated counts now on offer - all previous allocations have been cleared
                double[] allocCounts = sample.getPotentialAllocation(bucketGroup.getBucketRatios(), bucketGroup.getBucketIds(), bucketGroup.getBucketRatioRanges());
                double allocTotal = sumVector(allocCounts);

                if (allocTotal / sampleCount < reqAllocPercent)
                    continue;

                // add in descending order
                int index = 0;
                while(index < potentialAllocTotals.size())
                {
                    if(allocTotal > potentialAllocTotals.get(index))
                        break;

                    ++index;
                }

                bgIndexList.add(index, bgIndex);
                potentialAllocTotals.add(index, allocTotal);
                potentialAllocCounts.add(index, allocCounts);
            }

            if(bgIndexList.isEmpty())
                continue;

            int index = 0;

            int groupCount = bgIndexList.size();

            List<double[]> ratiosCollection = Lists.newArrayList();
            double[] potentialGroupContribs = new double[groupCount];

            final double[] sampleCounts = sample.getElevatedBucketCounts();
            final double[] countsNoise = sample.getCountRanges();
            int backgroundGroupIndex = -1;

            for (index = 0; index < groupCount; ++index)
            {
                int bgIndex = bgIndexList.get(index);
                final BucketGroup bucketGroup = mFinalBucketGroups.get(bgIndex);

                ratiosCollection.add(bucketGroup.getBucketRatios());

                potentialGroupContribs[index] = potentialAllocTotals.get(index);

                if(bucketGroup.equals(sample.getBackgroundGroup()))
                    backgroundGroupIndex = index;
            }

            sigOptim.initialise(sample.Id, sampleCounts, countsNoise, ratiosCollection, potentialGroupContribs);

            // each sample's background sig will remain in the list even if it drops below the required threshold
            sigOptim.setRequiredSig(backgroundGroupIndex);
            boolean validCalc = sigOptim.fitToSample(SAMPLE_ALLOCATED_PERCENT, MIN_GROUP_ALLOC_PERCENT_LOWER);

            if(!validCalc)
            {
                LOGGER.warn("sample({}) refit of {} sigs failed", sample.Id, groupCount);
                continue;
            }

            double[] newGroupContribs = new double[groupCount];
            copyVector(sigOptim.getContribs(), newGroupContribs);

            double fitAllocPerc = sigOptim.getAllocPerc();
            boolean useFittedAllocs = false;

            double fitAllocTotal = sumVector(newGroupContribs);

            if(fitAllocPerc > prevTotalAllocPerc)
            {
                LOGGER.debug(String.format("sample(%d) taking fitted allocation new(%.3f) vs prev(%.3f)",
                        sample.Id, fitAllocPerc, prevTotalAllocPerc));

                useFittedAllocs = true;
            }
            else
            {
                copyVector(potentialGroupContribs, newGroupContribs);
            }

            // finally allocate any group more than the min allocation percent
            List<Integer> sortedAllocIndices = getSortedVectorIndices(newGroupContribs, false);

            double maxAllocTotal = sumVector(potentialGroupContribs);

            for(index = 0; index < sortedAllocIndices.size(); ++index)
            {
                int itemIndex = sortedAllocIndices.get(index);
                Integer bgIndex = bgIndexList.get(itemIndex);
                final BucketGroup bucketGroup = mFinalBucketGroups.get(bgIndex);
                double potentialAlloc = potentialGroupContribs[itemIndex];
                double fitAlloc = newGroupContribs[itemIndex];

                if(fitAlloc/sampleCount < reqAllocPercent)
                    break;

                double[] allocCounts = potentialAllocCounts.get(itemIndex);

                if(useFittedAllocs)
                {
                    // override with the fitted contribution
                    final double[] bucketRatios = bucketGroup.getBucketRatios();
                    for(int b = 0; b < mBucketCount; ++b)
                    {
                        allocCounts[b] = fitAlloc * bucketRatios[b];
                    }
                }

                double prevAllocPerc = sample.getAllocPercent();
                double actualAlloc = sample.allocateBucketCounts(allocCounts, reqAllocPercent);
                double allocPerc = actualAlloc / sample.getElevatedCount();

                if(allocPerc >= reqAllocPercent)
                {
                    bucketGroup.addSample(sample.Id, allocCounts, false);
                    sample.addElevBucketGroup(bucketGroup, allocPerc);

                    LOGGER.debug(String.format("sample(%d) added to bg(%d) count(pot=%s fit=%s act=%s of %s) allocatedPerc(%.3f -> %.3f noise=%.3f) groupCount(%d)",
                            sample.Id, bucketGroup.getId(), sizeToStr(potentialAlloc), sizeToStr(fitAlloc), sizeToStr(actualAlloc), sizeToStr(sampleCount),
                            prevAllocPerc, sample.getAllocPercent(), sample.getAllocNoise()/sampleCount, sample.getElevBucketGroups().size()));
                }
            }

            String allocPercChange = "unch";

            if(sample.getAllocPercent() > prevTotalAllocPerc + 0.01)
                allocPercChange = "better";
            else if(sample.getAllocPercent() < prevTotalAllocPerc - 0.01)
                allocPercChange = "worse";

            LOGGER.debug(String.format("sample(%d) final fit: groups(%d prev=%d max=%d) %s allocation(prev=%.3f new=%.3f, act=%s of %s) proposed(%s) fit(%s)",
                    sample.Id, sample.getElevBucketGroups().size(), prevGroupCount, groupCount, allocPercChange, prevTotalAllocPerc, sample.getAllocPercent(),
                    sizeToStr(sample.getAllocatedCount()), sizeToStr(sample.getElevatedCount()), sizeToStr(maxAllocTotal), sizeToStr(fitAllocTotal)));
        }
    }

    private void calcBucketRatioRanges(final BucketGroup bucketGroup, List<Integer> sampleIds, List<Double> sampleAllocs)
    {
        List<Double> bucketRangePercent = Lists.newArrayList();

        if(sampleAllocs.size() != sampleIds.size())
            return;

        // calculate a percentage range for noise around each bucket
        List<Integer> groupBuckets = bucketGroup.getBucketIds();

        final double[] bucketRatios = bucketGroup.getBucketRatios();

        double totalAlloc = bucketGroup.getPotentialAllocation();

        for(Integer bucket : groupBuckets)
        {
            double percentTotal = 0;
            double ratioPercentTotal = 0;

            for(int samIndex = 0; samIndex < sampleIds.size(); ++samIndex)
            {
                Integer sampleId = sampleIds.get(samIndex);
                double proposedAlloc = sampleAllocs.get(samIndex);

                final SampleData sample = mSampleData.get(sampleId);
                double bucketRange = sample.calcBucketRange(bucket);
                double sampleFactor = proposedAlloc / totalAlloc;

                percentTotal += sampleFactor;
                ratioPercentTotal += sampleFactor * bucketRange;
            }

            double avgRange = capValue(ratioPercentTotal / percentTotal, 0, 0.5);

            bucketRangePercent.add(avgRange);

            LOGGER.debug(String.format("bg(%d) bucket(%d) ratio(%.4f) range(%.4f)",
                    bucketGroup.getId(), bucket, bucketRatios[bucket], avgRange));
        }

        bucketGroup.setBucketRatioRanges(bucketRangePercent);
    }

    private double minAllocPercent(final SampleData sample)
    {
        // require say 10% of the remaining unallocated count, but with an upper minimum limit
        return max(sample.getUnallocPercent() * MIN_GROUP_ALLOC_PERCENT, MIN_GROUP_ALLOC_PERCENT_LOWER);
    }

    private void purgeSimilarBucketGroups(final BucketGroup comparisonGroup)
    {
        final List<Integer> bucketIds = comparisonGroup.getBucketIds();

        int initGroupCount = mTopAllocBucketGroups.size();
        int bgIndex = 0;
        while(bgIndex < mTopAllocBucketGroups.size())
        {
            BucketGroup bucketGroup = mTopAllocBucketGroups.get(bgIndex);

            if(bucketGroup.hasAnyBucket(bucketIds))
            {
                mTopAllocBucketGroups.remove(bgIndex);
            }
            else
            {
                ++bgIndex;
            }
        }

        if(mTopAllocBucketGroups.size() < initGroupCount)
        {
            LOGGER.debug("reduced top bucket groups: {} -> {}", initGroupCount, mTopAllocBucketGroups.size());
        }
    }

    private void logWorstAllocatedSamples()
    {
        int maxWorstSamples = 100;
        List<int[]> worstAllocatedSamples = Lists.newArrayList();

        int SAMPLE_ID = 0;
        int UNALLOC_AMT = 1;

        for(final SampleData sample : mSampleData)
        {
            if(sample.isExcluded())
                continue;;

            if(sample.getAllocPercent() > 0.95)
                continue;

            int unallocTotal = (int)round((1 - sample.getAllocPercent()) * sample.getElevatedCount());
            int worstIndex = 0;
            for(; worstIndex < worstAllocatedSamples.size(); ++worstIndex)
            {
                if(unallocTotal > worstAllocatedSamples.get(worstIndex)[UNALLOC_AMT])
                    break;
            }

            if(worstIndex <= maxWorstSamples)
            {
                int[] worstData = { sample.Id, unallocTotal };
                worstAllocatedSamples.add(worstIndex, worstData);
            }
        }

        // log worst x samples
        for(int worstIndex = 0; worstIndex < min(worstAllocatedSamples.size(), maxWorstSamples); ++worstIndex)
        {
            final int[] worstData = worstAllocatedSamples.get(worstIndex);
            final SampleData sample = mSampleData.get(worstData[SAMPLE_ID]);

            LOGGER.debug(String.format("%d: worst sample(%d: %s) cancer(%s) unallocated(%s of %s, perc=%.3f) groupCount(%d)",
                    worstIndex, sample.Id, sample.getSampleName(), sample.getCancerType(), sizeToStr(worstData[UNALLOC_AMT]),
                    sizeToStr(sample.getElevatedCount()), (1 - sample.getAllocPercent()), sample.getElevBucketGroups().size()));
        }
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

        // clear out sample's group list since can be incorrect after merging - they will be added again below
        for (SampleData sample : mSampleData)
        {
            sample.clearElevBucketGroups();
        }

        for (final BucketGroup bucketGroup : mBucketGroups)
        {
            final List<Integer> sampleIds = bucketGroup.getSampleIds();

            if(sampleIds.isEmpty())
                continue;

            final List<Integer> bucketIds = bucketGroup.getBucketIds();
            List<Double> sampleCountTotals = Lists.newArrayList();

            // recalc the bucket count totals since merging did not attempt to handle double-counting sample's counts
            double[] bucketCounts = new double[mBucketCount];

            int lowProbCount = 0;
            double totalBmr = 0;

            for (Integer sampleId : sampleIds)
            {
                SampleData sample = mSampleData.get(sampleId);

                sample.addElevBucketGroup(bucketGroup);

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
                    LOGGER.debug(String.format("bg(%d) category(%s) count(%d of %d) perc(group=%.3f category=%.3f)",
                            bucketGroup.getId(), catName, catSampleCount, sampleCount, samplesPerc, catPerc));
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

    private void logBucketGroups(boolean verbose)
    {
        List<BucketGroup> bgList = mFinalBucketGroups;

        // sort by score before logging - for now log in order of creation (by highest allocation)
        // Collections.sort(bgList);

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

        double totalCount = mElevatedCount;

        for (int i = 0; i < maxToLog; ++i)
        {
            BucketGroup bucketGroup = bgList.get(i);

            if (bucketGroup.calcScore() < minScore)
                break;

            if (verbose)
            {
                String bucketIdsStr = "";

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

                final BucketFamily family = findBucketFamily(bucketGroup);

                double groupPerc = bucketGroup.getTotalCount() / totalCount;

                String linkData = "";
                if(mBucketGroupsLinked)
                {
                    linkData = String.format("family(%d) closestBg(%d css=%.4f) ",
                            family != null ? family.getId() : -1,
                            bucketGroup.getClosestBG() != null ? bucketGroup.getClosestBG().getId() : -1, bucketGroup.getClosestBGCss());
                }

                if(!bucketGroup.getTag().isEmpty())
                {
                    linkData += String.format("tag=%s ", bucketGroup.getTag());
                }

                LOGGER.debug(String.format("rank %d: bg(%d) %scancer(%s) score(%.0f purity=%.2f LF=%.1f) samples(%d) variants(avg=%s total=%s perc=%.3f) buckets(%d: %s) effects(%s)",
                        i, bucketGroup.getId(), linkData, bucketGroup.getCancerType(), bucketGroup.calcScore(), bucketGroup.getPurity(),
                        bucketGroup.getLoadFactor(), bucketGroup.getSampleIds().size(),
                        sizeToStr(bucketGroup.getAvgCount()), sizeToStr(bucketGroup.getTotalCount()), groupPerc,
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
            if(mBucketGroupFileWriter == null)
            {
                mBucketGroupFileWriter = getNewFile(mOutputDir, mOutputFileId + "_ba_group_data.csv");

                mBucketGroupFileWriter.write("Rank,BgId,Type,CancerType,Effects,SampleCount,BucketCount,MutLoad,PotentialAlloc,Purity");

                for(int i = 0; i < mBucketCount; ++i)
                {
                    mBucketGroupFileWriter.write(String.format(",%d", i));
                }

                mBucketGroupFileWriter.newLine();
            }

            BufferedWriter writer = mBucketGroupFileWriter;

            for (int bgIndex = 0; bgIndex < mFinalBucketGroups.size(); ++bgIndex)
            {
                BucketGroup bucketGroup = mFinalBucketGroups.get(bgIndex);

                final BucketFamily family = findBucketFamily(bucketGroup);

                String groupType = bucketGroup.getTag().equalsIgnoreCase("background") ? "BG" : "Elev";

                // writer.write(String.format("%s,%s", mAllocationRunId >= 0 ? mAllocationRunId : "Final", bucketGroup.isSelected() ? "TRUE" : "FALSE"));

                writer.write(String.format("%d,%d,%s,%s,%s,%d,%d,%.0f,%.0f,%.2f",
                        bgIndex, bucketGroup.getId(), groupType, bucketGroup.getCancerType(), bucketGroup.getEffects(),
                        bucketGroup.getSampleIds().size(), bucketGroup.getBucketIds().size(),
                        sumVector(bucketGroup.getBucketCounts()), bucketGroup.getPotentialAllocation(), bucketGroup.getPurity()));

                double[] bucketRatios = bucketGroup.getBucketRatios();

                for(int i = 0; i < mBucketCount; ++i)
                {
                    writer.write(String.format(",%.6f", bucketRatios[i]));
                }

                writer.newLine();
            }
        }
        catch(IOException exception)
        {
            LOGGER.error("failed to write output file: bucket groups");
        }
    }

    private void logOverallStats()
    {
        mAllocatedCount = 0;
        int fullyAllocated = 0;

        for(final SampleData sample : mSampleData)
        {
            mAllocatedCount += sample.getAllocatedCount();

            if(sample.getAllocPercent() >= SAMPLE_ALLOCATED_PERCENT)
                ++fullyAllocated;
        }

        LOGGER.debug(String.format("overall: samples(%d alloc=%d) groups(%d) counts: total(%s) background(%s perc=%.3f) elevated(%s perc=%.3f) allocated(%s perc=%.3f)",
                mSampleCount, fullyAllocated, mFinalBucketGroups.size(), sizeToStr(mTotalCount), sizeToStr(mBackgroundCount), mBackgroundCount/mTotalCount,
                sizeToStr(mElevatedCount), mElevatedCount/mTotalCount, sizeToStr(mAllocatedCount), mAllocatedCount/mElevatedCount));
    }

    private void logSampleResults()
    {
        final List<String> categories = mExtSampleData.getFieldNames();

        int fullyAllocCount = 0;
        int partiallyAllocCount = 0;
        int noMatchCount = 0;
        double countAllocated = 0; // the actual number of variants
        final double[][] elevData = mElevatedCounts.getData();

        for(int sampleId = 0; sampleId < mSampleCount; ++sampleId)
        {
            final SampleData sample = mSampleData.get(sampleId);

            if(sample.isExcluded())
                continue;;

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
                double noiseTotal = sumVector(sample.getCountRanges());
                double noiseAllocPerc = sample.getAllocNoise() / noiseTotal;

                LOGGER.debug(String.format("sample(%d: %s) %s allocated: groups(%d) buckets(%d unalloc=%d) count(total=%s bg=%s elev=%s alloc=%.3f noise=%.2f of %s) cancer(%s) effects(%d: %s)",
                        sampleId, sample.getSampleName(), fullyAllocated ? "fully" : "partially",
                        sample.getElevBucketGroups().size(), samBucketList.size(), sample.getUnallocBuckets().size(),
                        sizeToStr(sample.getTotalCount()), sizeToStr(bgTotal),
                        sizeToStr(sample.getElevatedCount()), sample.getAllocPercent(), noiseAllocPerc, sizeToStr(noiseTotal),
                        cancerType, effectsCount, effects));
                continue;
            }

            LOGGER.debug("sample({}: {}) unallocated with no match: buckets({}) count({} bg={} total={}) cancer({}) effects({}: {})",
                    sampleId, sample.getSampleName(), samBucketList.size(), sizeToStr(sample.getElevatedCount()), sizeToStr(bgTotal),
                    sizeToStr(sample.getTotalCount()), cancerType, effectsCount, effects);

            ++noMatchCount;
        }

        double percAllocated = countAllocated / mElevatedCount;

        LOGGER.debug(String.format("sample summary: total(%d) elev(%d) alloc(%d, %.3f of %s) partial(%d) unalloc(%d)",
                mSampleCount, mAllSampleBucketGroups.size(), fullyAllocCount, percAllocated, sizeToStr(mElevatedCount),
                partiallyAllocCount, noMatchCount));

        logOverallStats();
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
                    final BucketGroup bg1 = mFinalBucketGroups.get(mSigToBgMapping.get(sigId1));
                    final BucketGroup bg2 = mFinalBucketGroups.get(mSigToBgMapping.get(sigId2));

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
                final BucketGroup bucketGroup = mFinalBucketGroups.get(mSigToBgMapping.get(proposedSigId));

                LOGGER.debug(String.format("external ref sig(%d) matches bg(%d: ct=%s eff=%s samples=%d purity=%.2f) with css(%.4f)",
                        externalSigId, bucketGroup.getId(), bucketGroup.getCancerType(), bucketGroup.getEffects(),
                        bucketGroup.getSampleIds().size(), bucketGroup.getPurity(), result[CSSR_VAL]));
            }
        }
    }

    private void formBucketFamilies()
    {
        List<BucketGroup> bgList = mBucketGroups;

        Collections.sort(bgList);

        LOGGER.debug("creating families from {} elevated bucket groups", bgList.size());

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

            double familyPerc = bucketFamily.getTotalCount() / mElevatedCount;

            // final BucketGroup bucketGroup = mBucketGroups.get(mSigToBgMapping.get(proposedSigId));
            LOGGER.debug(String.format("family(%d) created from %d bucketGroups, buckets(%d) samples(%d) count(%s perc=%.3f) avgCss(%.4f) cancer(%s) effects(%s)",
                    bucketFamily.getId(), bucketFamily.getBucketGroups().size(), bucketFamily.getBucketIds().size(), bucketFamily.getSampleIds().size(),
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

            if(!mSpecificCancer.isEmpty() && !mSpecificCancer.equals(cancerType))
                continue;

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
        bucketGroup.setTag("background");
        mBackgroundGroups.add(bucketGroup);

        final List<Double> medianCounts = mBucketMediansMap.get(cancerType);

        if(medianCounts == null)
            return;

        double[] bucketRatios = listToArray(medianCounts);

        // for now add all samples to a single group per cancer type
        for (int index = 0; index < sampleIds.size(); ++index)
        {
            int sampleId = sampleIds.get(index);
            SampleData sample = mSampleData.get(sampleId);

            sample.setBackgroundGroup(bucketGroup);

            double[] sampleCounts = null;

            if(mNoBackgroundCounts)
            {
                sampleCounts = sample.getPotentialElevCounts(bucketRatios, fullBucketSet);
                double bgAllocTotal = sample.allocateBucketCounts(sampleCounts, 0);
                double bgAllocPerc = bgAllocTotal/sample.getTotalCount();
                sample.addElevBucketGroup(bucketGroup, bgAllocPerc);

                LOGGER.debug(String.format("sample(%d) added to background group(%d) alloc(%s perc=%.3f) ",
                        sampleId, bucketGroup.getId(), sizeToStr(bgAllocTotal), bgAllocPerc));
            }
            else
            {
                sampleCounts = mBackgroundCounts.getCol(sampleId);
            }

            bucketGroup.addSample(sampleId, sampleCounts, false);
        }

        // now counts are set, manually set the ratios (ie irrespective of the sample counts)
        bucketGroup.setBucketRatios(bucketRatios);
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

        for(int bgIndex = 0; bgIndex < mFinalBucketGroups.size(); ++bgIndex)
        {
            final BucketGroup bucketGroup = mFinalBucketGroups.get(bgIndex);

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

    private double[] extractBucketCountSubset(SampleData sample, final List<Integer> bucketSubset)
    {
        // extract the counts for the specified subset, leaving the rest zeroed
        double[] vec = new double[mBucketCount];

        final double[] unallocatedElevCounts = sample.getUnallocBucketCounts();

        for(Integer bucketId : bucketSubset)
        {
            vec[bucketId] = max(unallocatedElevCounts[bucketId], 0);
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

        for(final String cancerType : cancerTypes)
        {
            List<Integer> samplesList = mCancerSamplesMap.get(cancerType);

            if(samplesList.size() < minSamples)
            {
                minorsList.addAll(samplesList);
                mCancerSamplesMap.remove(cancerType);
            }
        }

        if(!minorsList.isEmpty())
        {
            mCancerSamplesMap.put(minors, minorsList);
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

            writer.write("SampleIndex,SampleId,CancerType,MutLoad,BackgroundCount,ElevCount,Status,AllocPerc,BgCount,BgIds,BgCancerType,BgEffects,ElevBuckets");

            // all unallocated bucket counts will be written out
            for(int i = 0; i < mBucketCount; ++i)
            {
                writer.write(String.format(",%d", i));
            }

            writer.newLine();

            for(final SampleData sample : mSampleData)
            {
                if(sample.isExcluded())
                    continue;;

                final List<Integer> sampleBuckets = sample.getElevatedBuckets();

                writer.write(String.format("%d,%s,%s", sample.Id, sample.getSampleName(), sample.getCancerType()));

                double bgTotal = sumVector(mBackgroundCounts.getCol(sample.Id));

                writer.write(String.format(",%.0f,%.0f,%.0f", sample.getTotalCount(), bgTotal, sample.getElevatedCount()));

                if(sampleBuckets == null)
                {
                    writer.write(",NoElevation,0,0,,,,");
                }
                else
                {
                    double allocPerc = sample.getAllocPercent();

                    String status = "Unalloc";
                    if (allocPerc >= 0.9)
                    {
                        status = "Alloc";
                    }
                    else if (allocPerc > 0.1)
                    {
                        status = "Partial";
                    }

                    final List<BucketGroup> bucketGroups = sample.getElevBucketGroups();
                    final List<Double> bgAllocPercents = sample.getGroupAllocPercents();

                    writer.write(String.format(",%s,%.3f,%d", status, allocPerc, bucketGroups.size()));

                    if (!bucketGroups.isEmpty())
                    {
                        String bgIdsStr = "";

                        for (int bgIndex = 0; bgIndex < bucketGroups.size(); ++bgIndex)
                        {
                            final BucketGroup bucketGroup = bucketGroups.get(bgIndex);

                            if (!bgIdsStr.isEmpty())
                                bgIdsStr += ";";

                            if(bgAllocPercents.size() == bucketGroups.size())
                                bgIdsStr += String.format("%d=%.3f", bucketGroup.getId(), bgAllocPercents.get(bgIndex));
                            else
                                bgIdsStr += bucketGroup.getId();
                        }

                        writer.write(String.format(",%s,%s,%s", bgIdsStr, bucketGroups.get(0).getCancerType(), bucketGroups.get(0).getEffects()));
                    }
                    else
                    {
                        writer.write(",,,");
                    }

                    String bucketIdsStr = sampleBuckets.toString().substring(1, sampleBuckets.toString().length() - 1); // remove []
                    if (sampleBuckets.size() > 1)
                        bucketIdsStr = bucketIdsStr.replaceAll(", ", ";");

                    writer.write(String.format(",%s", bucketIdsStr));
                }

                final double[] unallocCounts = sample.getUnallocBucketCounts();

                for(int i = 0; i < mBucketCount; ++i)
                {
                    writer.write(String.format(",%.0f", unallocCounts[i]));
                }

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


}
