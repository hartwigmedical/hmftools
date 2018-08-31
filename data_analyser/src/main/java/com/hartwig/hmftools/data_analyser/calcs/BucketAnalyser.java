package com.hartwig.hmftools.data_analyser.calcs;

import static java.lang.Math.abs;
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
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.convertToPercentages;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.copyVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.doubleToStr;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getDiffList;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getMatchingList;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getNewFile;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.listToArray;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sizeToStr;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sumVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sumVectors;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.vectorMultiply;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.writeMatrixData;
import static com.hartwig.hmftools.data_analyser.calcs.NmfConfig.NMF_REF_SIG_FILE;
import static com.hartwig.hmftools.data_analyser.types.GenericDataCollection.GD_TYPE_STRING;
import static com.hartwig.hmftools.data_analyser.types.NmfMatrix.redimension;
import static com.hartwig.hmftools.data_analyser.types.SampleData.PARTIAL_ALLOC_PERCENT;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.PerformanceCounter;
import com.hartwig.hmftools.data_analyser.loaders.GenericDataLoader;
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
    private NmfMatrix mPredefinedSigs;
    private boolean mFinalFitOnly; // using predefined sigs

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
    private List<BucketGroup> mSkippedBucketGroups;
    private List<Integer> mSkippedSamples;
    private List<Integer> mReassessSamples;
    private int mLastRunGroupCount;

    private boolean mNoBackgroundCounts; // whether to make a distinction between background and elevated counts
    private boolean mApplyNoise; // whether to factor Poisson noise into the sample counts and fits

    BufferedWriter mBucketGroupFileWriter;
    PerformanceCounter mPerfCounter;

    // config
    private String mOutputDir;
    private String mOutputFileId;
    private String mSpecificCancer;
    private int mSpecificSampleId; // purely for testing
    private double mHighCssThreshold; // CSS level for samples or groups to be consider similar
    private int mMaxProposedSigs;
    private int mMutationalLoadCap;
    private int mApplyPredefinedSigCount; // how many loaded sigs to apply prior to discovery (optimisation)
    private int mRunCount; //  number of iterations searching for potential bucket groups / sigs
    private int mRunId; // current run iteration
    private boolean mLogVerbose;
    private int mMinSampleAllocCount; // hard lower limit to allocate a sample to a group

    private static String BA_CSS_HIGH_THRESHOLD = "ba_css_high";
    private static String BA_MAX_PROPOSED_SIGS = "ba_max_proposed_sigs";
    private static String BA_CSS_SIG_THRESHOLD = "ba_css_proposed_sigs";
    private static String BA_SPECIFIC_CANCER = "ba_specific_cancer";
    private static String BA_RUN_COUNT = "ba_run_count";
    private static String BA_PREDEFINED_SIGS = "ba_predefined_sigs_file";
    private static String BA_PREDEFINED_SIG_APPLY_COUNT = "ba_predefined_sig_apply_count";
    private static String BA_MIN_SAM_ALLOC_COUNT = "ba_min_sample_alloc_count";
    private static String BA_MUT_LOAD_CAP = "ba_mut_load_cap";

    // constraint constants - consider moving any allocation-related ones to config to make them visible
    private static double SAMPLE_ALLOCATED_PERCENT = 0.995;
    private static double DOMINANT_CATEGORY_PERCENT = 0.7; /// mark a group with a category if X% of samples in it have this attribute (eg cancer type, UV)
    private static double MAX_ELEVATED_PROB = 1e-12;
    private int MIN_BUCKET_COUNT_OVERLAP = 3; // used for pairings of samples with reduced bucket overlap
    private static double MIN_DISCOVERY_SAMPLE_COUNT = 0.0001; // % of cohort to consider a pair of samples similar
    private static double PERMITTED_PROB_NOISE = 1e-4;
    private static double MIN_GROUP_ALLOC_PERCENT = 0.10; // only allocate a sample to a group if it takes at this much of the elevated count
    private static double MIN_GROUP_ALLOC_PERCENT_LOWER = 0.03; // hard lower limit
    private static double SKIP_ALLOC_FACTOR = 2; // skip adding a sample to a group if another candidate not chosen to allocate X times as much
    private static int MAX_CANDIDATE_GROUPS = 1500; // in place for speed and memory considerations
    private static double MAX_NOISE_TO_SAMPLE_RATIO = 5; // per sample, the total potential noise across all buckets cannot exceeds this multiple of variant total
    public static double MAX_NOISE_ALLOC_PERC = 0.5; // per sample, the max vs total variants which can be allocated to noise

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
    private boolean mHasErrors;

    public BucketAnalyser()
    {
        mOutputDir = "";
        mOutputFileId = "";
        mHasErrors = false;

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
        mSkippedSamples = Lists.newArrayList();
        mReassessSamples = Lists.newArrayList();
        mSkippedBucketGroups = Lists.newArrayList();
        mBucketGroupFileWriter = null;

        mHighCssThreshold = 0.995;
        mMaxProposedSigs = 0;
        mMutationalLoadCap = 0;
        mLogVerbose = false;
        mRunId = 0;
        mLastRunGroupCount = 0;
        mApplyPredefinedSigCount = 0;
        mFinalFitOnly = false;

        mPerfCounter = new PerformanceCounter("BucketAnalyser");

        mNoBackgroundCounts = true;
        mApplyNoise = true;

        mSampleWatchList = Lists.newArrayList();

        mSampleWatchList.add(1848);
        //mSampleWatchList.add(139);
        // mSampleWatchList.add(1013);
//        mSampleWatchList.add(2617);

        mSpecificSampleId = -1;
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(BA_EXT_SAMPLE_DATA_FILE, true, "Sample external data");
        options.addOption(BA_EXT_SAMPLE_DATA_FILE, true, "Sample external data");
        options.addOption(BA_CSS_HIGH_THRESHOLD, true, "Cosine sim for high-match test");
        options.addOption(BA_CSS_SIG_THRESHOLD, true, "Cosine sim for comparing proposed sigs");
        options.addOption(BA_MAX_PROPOSED_SIGS, true, "Maximum number of bucket groups to turn into proposed sigs");
        options.addOption(BA_SAMPLE_CALC_DATA_FILE, true, "Optional: file containing computed data per sample");
        options.addOption(BA_SPECIFIC_CANCER, true, "Optional: Only process this cancer type");
        options.addOption(BA_RUN_COUNT, true, "Number of search iterations");
        options.addOption(BA_PREDEFINED_SIGS, true, "Predefined sigs to use during discovery");
        options.addOption(BA_PREDEFINED_SIG_APPLY_COUNT, true, "How many predefined sigs to apply prior to discovery");
        options.addOption(BA_MIN_SAM_ALLOC_COUNT, true, "Min count to allocate a sample to a group");
        options.addOption(BA_MUT_LOAD_CAP, true, "Mutational load cap used in background counts calc");
    }

    public boolean initialise(GenericDataCollection collection, final CommandLine cmd)
    {
        mDataCollection = collection;
        mOutputFileId = cmd.getOptionValue(OUTPUT_FILE_ID);
        mOutputDir = cmd.getOptionValue(OUTPUT_DIR);
        mSpecificCancer = cmd.getOptionValue(BA_SPECIFIC_CANCER, "");

        mHighCssThreshold = Double.parseDouble(cmd.getOptionValue(BA_CSS_HIGH_THRESHOLD, "0.995"));
        mMaxProposedSigs = Integer.parseInt(cmd.getOptionValue(BA_MAX_PROPOSED_SIGS, "0"));
        mMutationalLoadCap = Integer.parseInt(cmd.getOptionValue(BA_MUT_LOAD_CAP, "10000"));
        mApplyPredefinedSigCount = Integer.parseInt(cmd.getOptionValue(BA_PREDEFINED_SIG_APPLY_COUNT, "0"));
        mMinSampleAllocCount = Integer.parseInt(cmd.getOptionValue(BA_MIN_SAM_ALLOC_COUNT, "1"));

        mRunCount = Integer.parseInt(cmd.getOptionValue(BA_RUN_COUNT, "25"));

        LOGGER.info(String.format("config: cssThreshold(%f) runCount(%d)",
                mHighCssThreshold, mRunCount));

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

        if(cmd.hasOption(BA_PREDEFINED_SIGS))
        {
            GenericDataCollection dataCollection = GenericDataLoader.loadFile(cmd.getOptionValue(BA_PREDEFINED_SIGS));
            mPredefinedSigs = DataUtils.createMatrixFromListData(dataCollection.getData());
            mPredefinedSigs.cacheTranspose();
            mFinalFitOnly = (mApplyPredefinedSigCount >= mPredefinedSigs.Cols);
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
            else
            {
                LOGGER.error("sampe({}) cannot find external data");
                mHasErrors = true;
            }

            if(!mSpecificCancer.isEmpty() && !sample.getCancerType().equals(mSpecificCancer))
            {
                sample.setExcluded(true);
            }
            else if(mSpecificSampleId >= 0 && sample.Id != mSpecificSampleId)
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

        return !mHasErrors;
    }

    public void run()
    {
        if(mHasErrors)
        {
            LOGGER.warn("failed to initialise, aborting run");
            return;
        }

        mPerfCounter. start();

        PerformanceCounter perfCounter = new PerformanceCounter("BucketMeanRatios");

        perfCounter.start("SplitCounts");
        calcBucketMedianData();
        calcSampleBackgroundCounts();
        splitSampleCounts();
        calcCountsNoise();
        collectElevatedSampleBuckets();

        // back-ground groups using expected counts
        formBackgroundBucketGroups();
        logOverallStats();
        perfCounter.stop();

        if(mApplyPredefinedSigCount > 0)
        {
            applyPredefinedSigs();
        }

        if(mFinalFitOnly)
            mRunCount = 0;

        for(mRunId = 0; mRunId < mRunCount; ++mRunId)
        {
            perfCounter.start(String.format("FindBucketGroup run %d", mRunId));

            clearBucketGroups();

            mLastRunGroupCount = mBucketGroups.size();

            LOGGER.debug("starting run({}) to find next top bucket group", mRunId);

            // search for precise bucket groups
            // formExactBucketGroups(false);

            // and then the best sub-groups
            formBucketGroupsFromSamplePairs(false);

            /*
            if(runId > 0)
            {
                formExactBucketGroups(true);
                formBucketGroupsFromSamplePairs(true);
            }
            */

            if(mBucketGroups.isEmpty())
            {
                LOGGER.debug("no new groups found for run({})", mRunId);
                break;
            }

            double prevAllocCount = mAllocatedCount;

            populateTopBucketGroups();

            // uncomment to log interim groups
            // analyseGroupsVsExtData(mTopAllocBucketGroups, false);
            // writeInterimBucketGroups();

            // allocated elevated counts to the best group
            BucketGroup nextBestGroup = allocateTopBucketGroup();

            if(nextBestGroup != null)
            {
                if(mFinalBucketGroups.contains(nextBestGroup))
                {
                    LOGGER.error("attempt to add bg({}) again", nextBestGroup.getId());
                    mHasErrors = true;
                    break;
                }

                logOverallStats();
                mFinalBucketGroups.add(nextBestGroup);
            }

            perfCounter.stop();

            analyseGroupsVsExtData(mFinalBucketGroups, false);
            logBucketGroups(true);

            if (nextBestGroup == null)
                break;

            if((mAllocatedCount - prevAllocCount) / mElevatedCount < 0.001) // eg 55K out of 55M
            {
                LOGGER.debug(String.format("breaking at run(%d) due to negligible allocPercChange(%s -> %s)",
                        mRunId, doubleToStr(prevAllocCount), doubleToStr(mAllocatedCount)));
                break;
            }

            if(mHasErrors)
            {
                LOGGER.warn("exiting run with errors");
                return;
            }
        }

        perfCounter.start("FinalFit");
        fitAllSamples();
        perfCounter.stop();

        perfCounter.start("AnalyseResults");
        analyseGroupsVsExtData(mFinalBucketGroups, true);
        logBucketGroups(true);
        logSampleResults();
        logWorstAllocatedSamples();
        logSimilarSampleContribs();
        perfCounter.stop();

        if(mMaxProposedSigs > 0)
        {
            createSignatures();
            compareSignatures();
            writeSampleContributions();
        }

        writeFinalBucketGroups();

        writeSampleData();
        // writeBackgroundSigs();
        // writeSampleMatrixData(mElevatedCounts, "_ba_elevated_counts.csv");
        writeSampleCalcData();

        perfCounter.logStats();

        finalise();

        mPerfCounter.stop();
        mPerfCounter.logStats();
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

        for(int s = 0; s < mSampleCount; ++s)
        {
            mSampleBgAllocations.add(s, 0.0);
        }

        final double[][] scData = mSampleCounts.getData();

        for (Map.Entry<String, List<Double>> entry : mBucketMediansMap.entrySet())
        {
            final String type = entry.getKey();

            List<Double> medianRatios = entry.getValue();

            List<Integer> sampleIds = mCancerSamplesMap.get(type);

            final double[] bucketRatios = listToArray(medianRatios);

            for (Integer sampleId : sampleIds)
            {
                double[] sampleBgCounts = new double[mBucketCount];

                for (int i = 0; i < mBucketCount; ++i)
                {
                    int expectedCount = (int) round(bucketRatios[i] * min(mSampleTotals[sampleId], mMutationalLoadCap));
                    sampleBgCounts[i] = min(expectedCount, scData[i][sampleId]);
                }

                double optimalBgCount = calcBestFitWithinProbability(sampleId, bucketRatios, sampleBgCounts, 0.99, 0.01);
                mSampleBgAllocations.set(sampleId, optimalBgCount);
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
                if(!mNoBackgroundCounts)
                {
                    int backgroundCount = (int) bgData[i][j];

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
                }

                int elevatedCount = (int)elevData[i][j];

                // compute a range for Poisson noise around this elevated count
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

        mPermittedElevRange.cacheTranspose();
        mPermittedBgRange.cacheTranspose();

        /*
        // TEMP to write out sample noise values
        try
        {
            BufferedWriter writer = getNewFile(mOutputDir, mOutputFileId + "_ba_sample_noise.csv");

            List<String> sampleNames = mDataCollection.getFieldNames();

            int i = 0;
            for(; i < sampleNames.size()-1; ++i)
            {
                writer.write(String.format("%s,", sampleNames.get(i)));
            }
            writer.write(String.format("%s", sampleNames.get(i)));

            writer.newLine();

            writeMatrixData(writer, mPermittedElevRange, true);

            writer.close();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing sample noise file");
        }
        */
    }

    private void collectElevatedSampleBuckets()
    {
        mAllSampleBucketGroups = new HashMap();

        int totalCount = 0;
        double[][] probData = mBucketProbs.getData();

        double[] countRanges = new double[mBucketCount];

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

            if(mApplyNoise)
            {
                copyVector(mPermittedElevRange.getCol(i), countRanges);

                // ensure that noise counts aren't disproportionately large compared with the actual counts
                double noiseCountsTotal = sumVector(countRanges);
                double sampleNoiseTotal = CosineSim.calcPoissonRangeGivenProb((int)round(mSampleTotals[i]), PERMITTED_PROB_NOISE);

                if(noiseCountsTotal > MAX_NOISE_TO_SAMPLE_RATIO * sampleNoiseTotal)
                {
                    double reduceRatio = sampleNoiseTotal * MAX_NOISE_TO_SAMPLE_RATIO / noiseCountsTotal;

                    for (int j = 0; j < mBucketCount; ++j)
                    {
                        countRanges[j] = round(countRanges[j] * reduceRatio);
                    }
                }
            }

            sample.setElevatedBucketCounts(mElevatedCounts.getCol(i), countRanges);

            if (!bucketList.isEmpty())
            {
                sample.setElevatedBuckets(bucketList);
                mAllSampleBucketGroups.put(i, bucketList);
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

        BucketGroup bucketGroup = new BucketGroup(mNextBucketId++);
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

            if(mSampleWatchList.contains(sampleId))
            {
                // LOGGER.debug("spec sample");
            }

            double[] sampleCounts = null;

            if(mNoBackgroundCounts)
            {
                sampleCounts = sample.getPotentialElevCounts(bucketRatios, fullBucketSet);
                double bgAllocTotal = sample.allocateBucketCounts(sampleCounts, 0);
                sample.addElevBucketGroup(bucketGroup, sample.getAllocPercent());

                LOGGER.debug(String.format("sample(%d) added to background group(%d) alloc(%s perc=%.3f of %s) noise(%s, %.3f of %s)",
                        sampleId, bucketGroup.getId(), sizeToStr(bgAllocTotal), sample.getAllocPercent(), sizeToStr(sample.getTotalCount()),
                        sizeToStr(sample.getAllocNoise()), sample.getNoisePerc(), sizeToStr(sample.getNoiseTotal())));
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

    private void formBucketGroupsFromSamplePairs(boolean usePartialUnallocated)
    {
        // create unique groups from any subset of 2 samples' buckets if they are a close match
        // samples aren't necessarily allocated, only groups are made for the time-being
        LOGGER.debug("forming bucket groups from sample-pair reduced buckets");

        String bgTag = usePartialUnallocated ? "subset-partial" : "subset";

        double[] sc1 = new double[mBucketCount];
        double[] sc2 = new double[mBucketCount];

        double minSampleCount = MIN_DISCOVERY_SAMPLE_COUNT * mTotalCount;

        int groupsCreated = 0;

        for (int samIndex1 = 0; samIndex1 < mSampleCount; ++samIndex1)
        {
            SampleData sample1 = mSampleData.get(samIndex1);

            if(sample1.isExcluded())
                continue;

            if(!mReassessSamples.isEmpty() && !mReassessSamples.contains(samIndex1))
                continue;

            if(usePartialUnallocated && sample1.getAllocPercent() < PARTIAL_ALLOC_PERCENT)
                continue;

            if(mSampleWatchList.contains(samIndex1))
            {
                //LOGGER.debug("spec sample");
            }

            double reqSam1AllocPercent = minAllocPercent(sample1, false);

            if(sample1.getUnallocPercent() < reqSam1AllocPercent)
                continue;

            if(sample1.getTotalCount() < minSampleCount)
                continue;

            final List<Integer> bl1 = usePartialUnallocated ? sample1.getElevatedBuckets() : sample1.getUnallocBuckets();

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

                if(mSampleWatchList.contains(samIndex1) && mSampleWatchList.contains(samIndex2))
                {
                    //LOGGER.debug("spec sample");
                }

                if(usePartialUnallocated && sample2.getAllocPercent() < PARTIAL_ALLOC_PERCENT)
                    continue;

                double reqSam2AllocPercent = minAllocPercent(sample2, false);

                if(sample2.getUnallocPercent() < reqSam2AllocPercent)
                    continue;

                if(sample2.getTotalCount() < minSampleCount)
                    continue;

                final List<Integer> bl2 = usePartialUnallocated ? sample2.getElevatedBuckets() : sample2.getUnallocBuckets();

                if (bl2.isEmpty())
                    continue;

                List<Integer> commonBuckets = getMatchingList(bl1, bl2);
                int commonBucketCount = commonBuckets.size();

                if (commonBucketCount < 2)
                    continue;

                sample1.populateBucketCountSubset(sc1, commonBuckets, usePartialUnallocated);
                sample2.populateBucketCountSubset(sc2, commonBuckets, usePartialUnallocated);

                final double[] sam1ElevCounts = usePartialUnallocated ? sample1.getPartialUnallocBucketCounts() : sample1.getUnallocBucketCounts();
                final double[] sam2ElevCounts = usePartialUnallocated ? sample2.getPartialUnallocBucketCounts() : sample2.getUnallocBucketCounts();

                double sam1ElevTotal = sumVector(sc1);
                double sam2ElevTotal = sumVector(sc2);
                double elevatedTotal = sam1ElevTotal + sam2ElevTotal;

                if(sam1ElevTotal/sample1.getElevatedCount() < reqSam1AllocPercent || sam2ElevTotal/sample2.getElevatedCount() < reqSam2AllocPercent)
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
            BucketGroup bucketGroup = new BucketGroup(mNextBucketId++);
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

            if(mLogVerbose)
            {
                LOGGER.debug(String.format("added bg(%d) samples(%d and %d) with buckets(%d) css(%.4f) allocCalcTotal(%s)",
                        bucketGroup.getId(), samIndex1, maxOtherSample, maxSharedBuckets.size(), maxCss, sizeToStr(maxAllocaTotal)));
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

        LOGGER.debug("created {} sample-pair {} groups", groupsCreated, bgTag);
    }

    private void populateTopBucketGroups()
    {
        mTopAllocBucketGroups.clear();

        LOGGER.debug("finding top potential bucket group from count({} new={})", mBucketGroups.size(), mBucketGroups.size() - mLastRunGroupCount);

        int maxCandidateGroups = MAX_CANDIDATE_GROUPS;

        SigContributionOptimiser sigOptim = new SigContributionOptimiser(mBucketCount, false, SAMPLE_ALLOCATED_PERCENT, true);

        int exceededOnSoloAlloc = 0;
        int exceededOnUnalloc = 0;
        int exceededOnFit = 0;
        int skippedRetry = 0;

        boolean useLastRunLogic = mLastRunGroupCount > 0;

        // first clear all existing allocations of samples to groups and vice versa
        for (int bgIndex = 0; bgIndex < mBucketGroups.size(); ++bgIndex)
        {
            BucketGroup bucketGroup = mBucketGroups.get(bgIndex);

            if(!useLastRunLogic)
            {
                bucketGroup.clearSamples();
                bucketGroup.resetPotentialAllocation();
            }

            double[] bgRatios = new double[mBucketCount];
            copyVector(bucketGroup.getBucketRatios(), bgRatios); // won't be recomputed as sample counts are added

            final List<Integer> groupBuckets = bucketGroup.getBucketIds();

            for (int sampleId = 0; sampleId < mSampleCount; ++sampleId)
            {
                final SampleData sample = mSampleData.get(sampleId);

                if(sample.isExcluded())
                    continue;

                if(mSampleWatchList.contains(sampleId) && mSampleWatchList.contains(bucketGroup.getId())) //
                {
                    // LOGGER.debug("spec sample");
                }

                double reqAllocPercent = minAllocPercent(sample, false);
                boolean exceedsMinAllocPerc = false;
                double[] allocCounts = null;
                double allocCountTotal = 0;
                double allocPercent = 0;

                if(useLastRunLogic)
                {
                    // pre-existing bucket groups (ie those not just proposed) and samples just not allocated can be left alone
                    if (bgIndex < mLastRunGroupCount && !mReassessSamples.contains(sampleId))
                    {
                        // look for an existing allocation in this group
                        if (bucketGroup.hasSample(sampleId))
                        {
                            allocCountTotal = bucketGroup.getSampleCount(sampleId);

                            if(allocCountTotal/sample.getElevatedCount() < reqAllocPercent)
                            {
                                LOGGER.error(String.format("sample(%d) part of existing bg(%d) with alloc(%s perc=%.3f) too low",
                                        sampleId, bucketGroup.getId(), sizeToStr(allocCountTotal), allocCountTotal/sample.getElevatedCount()));
                                mHasErrors = true;
                            }
                        }
                        else
                        {
                            // no point trying again
                        }

                        ++skippedRetry;
                        continue;
                    }
                }

                final List<Integer> samBuckets = sample.getElevatedBuckets();

                if (samBuckets.isEmpty() || sample.getUnallocPercent() < reqAllocPercent)
                    continue;

                // optimisation: check whether the buckets for this group and sample
                // could possibly exceed the min % threshold with a perfect fit, otherwise skip it
                final double[] bgSampleElevCounts = sample.getPotentialElevCounts(bgRatios, groupBuckets);
                double maxPotentialPerc = sumVector(bgSampleElevCounts) / sample.getElevatedCount();

                if(maxPotentialPerc < reqAllocPercent)
                    continue;

                // if the potential allocation taking no existing allocations into account would increase
                // a sample's overall allocation by more than the upper threshold, it must satisfy the test
                // to be added to this candidate
                if(maxPotentialPerc - sample.getAllocPercent() >= reqAllocPercent)
                {
                    exceedsMinAllocPerc = true;
                    ++exceededOnSoloAlloc;

                    allocPercent = maxPotentialPerc - sample.getAllocPercent();
                    allocCounts = bgSampleElevCounts;
                    allocCountTotal = sumVector(allocCounts);
                }
                else
                {
                    allocCounts = sample.getPotentialUnallocCounts(bgRatios, groupBuckets);
                    allocCountTotal = sumVector(allocCounts);
                    allocPercent = allocCountTotal / sample.getElevatedCount();

                    exceedsMinAllocPerc = allocPercent >= reqAllocPercent;

                    if(exceedsMinAllocPerc)
                        ++exceededOnUnalloc;
                }

                if (!exceedsMinAllocPerc)
                {
                    if(sample.getElevBucketGroups().isEmpty())
                        continue;

                    boolean sampleLog = false;
//                    if(mSampleWatchList.contains(sampleId) && mSampleWatchList.contains(bucketGroup.getId()) && sample.getElevBucketGroups().size() == 2)
//                    {
//                        // LOGGER.debug("spec sample");
//                        sampleLog = true;
//                    }

                    // see if a fit with sig along with all the other allocated one for this sample would then meet the min % threshold
                    // it is the overall change to the sample's allocation that is tested, not just this proposed group's contribution
                    List<double[]> ratiosCollection = Lists.newArrayList();

                    // intentionally add the candidate group before the others to give it best chance of allocation
                    // until the optim routine can do this better
                    // int candidateSigIndex = 0;

                    for (final BucketGroup samGroup : sample.getElevBucketGroups())
                    {
                        ratiosCollection.add(samGroup.getBucketRatios());
                    }

                    ratiosCollection.add(bgRatios);

                    double[] prevContribs = new double[ratiosCollection.size()];
                    int candidateSigIndex = prevContribs.length - 1;

                    sigOptim.initialise(sample.Id, sample.getElevatedBucketCounts(), sample.getCountRanges(), ratiosCollection, prevContribs, reqAllocPercent, mMinSampleAllocCount);
                    sigOptim.setLogVerbose(sampleLog);
                    sigOptim.setTargetSig(candidateSigIndex);
                    boolean validCalc = sigOptim.fitToSample();

                    if (!validCalc) // couldn't reach the required percent for this canidate sig
                    {
                        LOGGER.warn("sample({}) fit with existing sigs failed", sample.Id);
                        mHasErrors = true;
                        continue;
                    }

                    double fitAllocTotal = sigOptim.getContribTotal();
                    allocCountTotal = fitAllocTotal - (sample.getAllocatedCount() + sample.getAllocNoise());
                    allocPercent = allocCountTotal / sample.getElevatedCount();

                    if (allocPercent < reqAllocPercent || allocCountTotal < mMinSampleAllocCount)
                        continue;

                    // turn the fitted contribution into new counts
                    double candidateAlloc = sigOptim.getContribs()[candidateSigIndex];

                    if(candidateAlloc < allocCountTotal * 0.99)
                        continue; // shouldn't happen but in case a reshuffle of sigs didn't actually include the new sig

                    for(int b = 0; b < mBucketCount; ++b)
                    {
                        allocCounts[b] = bgRatios[b] * candidateAlloc;
                    }

                    exceedsMinAllocPerc = true;
                    ++exceededOnFit;
                }

                bucketGroup.addPotentialAllocation(allocCountTotal);
                bucketGroup.addPotentialAdjAllocation(allocCountTotal * allocPercent);
                bucketGroup.addSample(sampleId, allocCounts); // causes a recalc of ratios
            }

            if(bucketGroup.getPotentialAllocation() == 0)
                continue;

            // store in order since the top one(s) will be allocated first
            int findBgIndex = 0;
            while(findBgIndex < mTopAllocBucketGroups.size())
            {
                if(bucketGroup.getPotentialAdjAllocation() >= mTopAllocBucketGroups.get(findBgIndex).getPotentialAdjAllocation())
                    break;

                ++findBgIndex;
            }

            if(bgIndex < maxCandidateGroups || maxCandidateGroups == 0)
            {
                mTopAllocBucketGroups.add(findBgIndex, bucketGroup);
            }

            if(maxCandidateGroups > 0 && mTopAllocBucketGroups.size() > maxCandidateGroups)
            {
                mTopAllocBucketGroups.remove(mTopAllocBucketGroups.size() - 1);
            }
        }

        if(mTopAllocBucketGroups.isEmpty())
        {
            LOGGER.debug("no max potential group found");
        }

        LOGGER.debug("found {} top bucket groups, method(solo={} unalloc={} fit={} skipped={})",
                mTopAllocBucketGroups.size(), exceededOnSoloAlloc, exceededOnUnalloc, exceededOnFit, skippedRetry);
    }

    private BucketGroup allocateTopBucketGroup()
    {
        if (mTopAllocBucketGroups.isEmpty())
            return null;

        mReassessSamples.clear();

        BucketGroup topBucketGroup = mTopAllocBucketGroups.get(0);

        LOGGER.debug(String.format("top bg(%d) with buckets(%d) samples(%d) potential allocation(%s adj=%s)",
                topBucketGroup.getId(), topBucketGroup.getBucketIds().size(), topBucketGroup.getSampleIds().size(),
                sizeToStr(topBucketGroup.getPotentialAllocation()), sizeToStr(topBucketGroup.getPotentialAdjAllocation())));

        List<Integer> topSampleIds = Lists.newArrayList();
        topSampleIds.addAll(topBucketGroup.getSampleIds());

        List<Double> sampleAllocTotals = Lists.newArrayList();
        sampleAllocTotals.addAll(topBucketGroup.getSampleCountTotals());

        List<double[]> sampleCounts = Lists.newArrayList();
        sampleCounts.addAll(topBucketGroup.getSampleCounts());

        refineTopBucketGroup(topBucketGroup, topSampleIds, sampleAllocTotals, sampleCounts);

        final List<Integer> groupBuckets = topBucketGroup.getBucketIds();
        // final List<Double> bucketRatioRanges = topBucketGroup.getBucketRatioRanges();

        // now allocate samples to this top group
        topBucketGroup.clearSamples();

        double totalAlloc = 0;

        List<Integer> skippedSamples = Lists.newArrayList();

        for(int samIndex = 0; samIndex < topSampleIds.size(); ++samIndex)
        {
            Integer sampleId = topSampleIds.get(samIndex);
            double newAllocTotal = sampleAllocTotals.get(samIndex);
            double[] sampleCountAllocations = sampleCounts.get(samIndex);

            final SampleData sample = mSampleData.get(sampleId);

            if(mSampleWatchList.contains(sampleId))
            {
                //LOGGER.debug("spec sample");
            }

            BucketGroup maxOtherGroup = getSampleMaxAllocationGroup(sample, topBucketGroup);
            double maxOtherGroupAlloc = maxOtherGroup != null ? maxOtherGroup.getSampleCount(sampleId) : 0;
            double maxFinalGroupAlloc = getMaxAllocationAgainstFinalGroups(sample);
            double maxOtherAlloc = max(maxOtherGroupAlloc, maxFinalGroupAlloc);

            if(maxOtherAlloc > SKIP_ALLOC_FACTOR * newAllocTotal && maxOtherAlloc >= mMinSampleAllocCount)
            {
                if(maxOtherGroupAlloc > maxFinalGroupAlloc)
                {
                    LOGGER.debug(String.format("sample(%d) skipped bg(%d) alloc(%s) for other candidate bg(%d alloc=%s)",
                            sampleId, topBucketGroup.getId(), sizeToStr(newAllocTotal), maxOtherGroup.getId(), sizeToStr(maxOtherGroupAlloc)));

                    if (!mSkippedBucketGroups.contains(maxOtherGroup))
                    {
                        mSkippedBucketGroups.add(maxOtherGroup);
                    }
                }
                else
                {
                    LOGGER.debug(String.format("sample(%d) skipped bg(%d) alloc(%s) for final bg alloc(%s)",
                            sampleId, topBucketGroup.getId(), sizeToStr(newAllocTotal), sizeToStr(maxFinalGroupAlloc)));
                }

                skippedSamples.add(sampleId);
                continue;
            }

            if(newAllocTotal < mMinSampleAllocCount)
                continue;

            double reqAllocPercent = minAllocPercent(sample, false);

            // apply to the remaining unallocated elevated counts for this sample
            double actualAlloc = sample.allocateBucketCounts(sampleCountAllocations, reqAllocPercent);
            double allocPerc = actualAlloc / sample.getElevatedCount();

            if(allocPerc >= reqAllocPercent)
            {
                topBucketGroup.addSample(sampleId, sampleCountAllocations, false);
                sample.addElevBucketGroup(topBucketGroup, allocPerc);
                totalAlloc += actualAlloc;

                LOGGER.debug(String.format("sample(%d) added to bg(%d) buckets(grp=%d sam=%d unalloc=%d) count(proposed=%s act=%s of %s) allocatedPerc(+%.3f -> %.3f) noise(%s %.3f/%.3f) groupCount(%d)",
                        sampleId, topBucketGroup.getId(), groupBuckets.size(), sample.getElevatedBuckets().size(), sample.getUnallocBuckets().size(),
                        sizeToStr(newAllocTotal), sizeToStr(actualAlloc), sizeToStr(sample.getElevatedCount()), sample.lastAllocPercChange(), sample.getAllocPercent(),
                        sizeToStr(sample.getAllocNoise()), sample.getNoisePerc(), sample.getNoiseOfTotal(), sample.getElevBucketGroups().size()));
            }
        }

        LOGGER.debug(String.format("new top bg(%d) added %d samples, totalAllocatedCount(%s)",
                topBucketGroup.getId(), topBucketGroup.getSampleIds().size(), sizeToStr(totalAlloc)));

        // recheck any sample previously skipped against the fixed groups so see if they can now be allocated
        checkSkippedSamples(skippedSamples, topBucketGroup.getSampleIds());

        // for now clear all top groups to force a reassessment, since after allocation of elevated counts
        // to the top bucket, the similarities could increase across more buckets and samples
        mTopAllocBucketGroups.clear();

        if(topBucketGroup.getSampleIds().isEmpty())
            return null;

        mReassessSamples.addAll(topBucketGroup.getSampleIds());

        if(mSkippedBucketGroups.contains(topBucketGroup))
        {
            LOGGER.debug("previously skipped group({}) now added to final list", topBucketGroup.getId());
            mSkippedBucketGroups.remove(topBucketGroup);
        }

        return topBucketGroup;
    }

    private void refineTopBucketGroup(BucketGroup topBucketGroup, List<Integer> topSampleIds, List<Double> sampleAllocTotals, List<double[]> sampleCounts)
    {
        // walk down list, recalculating the ratios, expanding their range and comparing to groups further down
        final double[] topBucketRatios = topBucketGroup.getBucketRatios();
        double[] avgBucketRatios = new double[mBucketCount];
        copyVector(topBucketRatios, avgBucketRatios);

        int maxComparisons = 100;
        double topPotentialAlloc = topBucketGroup.getPotentialAllocation();
        List<Integer> candidateNewBuckets = Lists.newArrayList();
        List<BucketGroup> similarGroups = Lists.newArrayList();

        for (int bgIndex = 1; bgIndex < min(mTopAllocBucketGroups.size(), maxComparisons); ++bgIndex)
        {
            BucketGroup bucketGroup = mTopAllocBucketGroups.get(bgIndex);

            double[] bucketRatios = bucketGroup.getBucketRatios();

            double groupCss = calcCSS(topBucketRatios, bucketRatios);
            if (groupCss < mHighCssThreshold)
                continue;

            List<Integer> deficientBuckets = getDiffList(topBucketGroup.getBucketIds(), bucketGroup.getBucketIds());

            if(!deficientBuckets.isEmpty())
            {
                similarGroups.add(bucketGroup);
                continue;
            }

            // don't merge the group if there are no new samples to add
            List<Integer> newSamples = getDiffList(bucketGroup.getSampleIds(), topSampleIds);

            if(newSamples.isEmpty())
            {
                similarGroups.add(bucketGroup);
                continue;
            }

            double groupToTopRatio = bucketGroup.getPotentialAllocation() / topPotentialAlloc;
            vectorMultiply(bucketRatios, groupToTopRatio);

            for (Integer bucket : topBucketGroup.getBucketIds())
            {
                avgBucketRatios[bucket] += bucketRatios[bucket];
            }

            List<Integer> extraBuckets = getDiffList(bucketGroup.getBucketIds(), topBucketGroup.getBucketIds());

            for (Integer bucket : extraBuckets)
            {
                if (!candidateNewBuckets.contains(bucket))
                    candidateNewBuckets.add(bucket);
            }

            // merge in any difference in samples
            final List<Integer> bgSamples = bucketGroup.getSampleIds();
            final List<Double> bgSampleTotals = bucketGroup.getSampleCountTotals();
            final List<double[]> bgSampleCounts = bucketGroup.getSampleCounts();

            int samplesAdded = 0;
            for (int samIndex = 0; samIndex < bgSamples.size(); ++samIndex)
            {
                Integer sampleId = bgSamples.get(samIndex);

                if (newSamples.contains(sampleId))
                {
                    topSampleIds.add(sampleId);
                    sampleAllocTotals.add(bgSampleTotals.get(samIndex));
                    sampleCounts.add(bgSampleCounts.get(samIndex));
                    ++samplesAdded;
                }
            }

            LOGGER.debug(String.format("top bg(%d) merged with bg(%d) alloc(%s adj=%s) css(%.4f) samples(bg1=%d bg2=%d added=%d) diffBuckets(%d)",
                    topBucketGroup.getId(), bucketGroup.getId(), sizeToStr(bucketGroup.getPotentialAllocation()),
                    sizeToStr(bucketGroup.getPotentialAdjAllocation()), groupCss, topSampleIds.size(), bgSamples.size(), samplesAdded, extraBuckets.size()));

            // remove this from the candidate set so it doesn't impact the skipped-sample logic
            similarGroups.add(bucketGroup);
        }

        for (final BucketGroup candidateBg : mBucketGroups)
        {
            if (similarGroups.contains(candidateBg))
            {
                mBucketGroups.remove(candidateBg);
                break;
            }
        }

        // test out any candidate additional buckets
        // testCandidateExtraBuckets(topBucketGroup, avgBucketRatios, candiateAdditionalBuckets, topSampleIds, sampleAllocTotals)

        // convert back to percentages
        convertToPercentages(avgBucketRatios);

        List<SampleData> samples = Lists.newArrayList();
        List<double[]> proposedAllocs = Lists.newArrayList();
        for(int samIndex = 0; samIndex < topSampleIds.size(); ++samIndex)
        {
            int sampleId = topSampleIds.get(samIndex);
            samples.add(mSampleData.get(sampleId));
            proposedAllocs.add(sampleCounts.get(samIndex));
        }

        SigOptimiser sigOptim = new SigOptimiser(topBucketGroup.getId(), samples, proposedAllocs, avgBucketRatios, candidateNewBuckets);
        sigOptim.setLogVerbose(true);
        boolean validCalc = sigOptim.optimiseBucketRatios();

        if(validCalc && sigOptim.hasChanged())
        {
            copyVector(sigOptim.getFittedRatios(), avgBucketRatios);

            final List<Integer> newBuckets = sigOptim.getNewBuckets();
            for (Integer newBucket : newBuckets)
            {
                topBucketGroup.addBucket(newBucket, false);
            }

            // take the new samples allocations as well
            for (int samIndex = 0; samIndex < topSampleIds.size(); ++samIndex)
            {
                sampleAllocTotals.set(samIndex, sigOptim.getRevisedSampleAlloc(samIndex));
                sampleCounts.set(samIndex, sigOptim.getRevisedSampleAllocCounts(samIndex));
            }
        }

        topBucketGroup.setBucketRatios(avgBucketRatios);
    }

    private BucketGroup getSampleMaxAllocationGroup(final SampleData sample, final BucketGroup excludeGroup)
    {
        double maxAlloc = 0;

        // test against the candidate list of bucket groups
        BucketGroup bestGroup = null;

        for (final BucketGroup bucketGroup : mBucketGroups)
        {
            if (bucketGroup == excludeGroup)
                continue;

            double alloc = bucketGroup.getSampleCount(sample.Id);
            if (alloc > maxAlloc)
            {
                maxAlloc = alloc;
                bestGroup = bucketGroup;
            }
        }

        return bestGroup;
    }

    private double getMaxAllocationAgainstFinalGroups(final SampleData sample)
    {
        double maxAlloc = 0;

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

        return maxAlloc;
    }

    private void checkSkippedSamples(List<Integer> newSkippedSamples, List<Integer> newAddedSamples)
    {
        // for samples previously skipped, check if can now be allocated to one or more of the final bucket groups
        //if(mFinalBucketGroups.isEmpty())
        //    return;

        if(newSkippedSamples.isEmpty() && mSkippedSamples.isEmpty())
            return;

        int initSkippedSamples = mSkippedSamples.size();
        int allocatedSamples = 0;
        int reskippedSamples = 0;

        // first remove the ones just added to a group
        for(Integer sampleId : newAddedSamples)
        {
            if(mSkippedSamples.contains(sampleId))
                mSkippedSamples.remove(sampleId);
        }

        // now check whether any others previously skipped could now be added to one of the final groups
        if(!mFinalBucketGroups.isEmpty())
        {
            int samIndex = 0;
            while (samIndex < mSkippedSamples.size())
            {
                Integer sampleId = mSkippedSamples.get(samIndex);

                if (newSkippedSamples.contains(sampleId))
                {
                    ++samIndex;
                    continue;
                }

                SampleData sample = mSampleData.get(sampleId);

                // first check against the candidate groups
                double maxOtherGroupAlloc = 0;
                double maxOtherGroupAllocPerc = 0;
                BucketGroup maxOtherGroup = getSampleMaxAllocationGroup(sample, null);

                if (maxOtherGroup != null)
                {
                    maxOtherGroupAlloc = maxOtherGroup.getSampleCount(sample.Id);
                    maxOtherGroupAllocPerc = maxOtherGroupAlloc / sample.getElevatedCount();
                }

                double reqAllocPercent = minAllocPercent(sample, false);

                boolean wasAllocated = false;

                // now also check against the final groups
                for (final BucketGroup bucketGroup : mFinalBucketGroups)
                {
                    if (bucketGroup.hasSample(sampleId))
                        continue;

                    double[] allocCounts = sample.getPotentialUnallocCounts(bucketGroup.getBucketRatios(), bucketGroup.getBucketIds());

                    double proposedAllocTotal = sumVector(allocCounts);
                    double proposedAllocPerc = proposedAllocTotal / sample.getElevatedCount();

                    if (proposedAllocPerc < reqAllocPercent || proposedAllocTotal < mMinSampleAllocCount)
                        continue;

                    // continue to hold out if there is a better candidate group yet to come
                    if (maxOtherGroupAlloc > SKIP_ALLOC_FACTOR * proposedAllocTotal)
                    {
                        ++reskippedSamples;
                        continue;
                    }

                    // apply to the remaining unallocated elevated counts for this sample
                    double actualAlloc = sample.allocateBucketCounts(allocCounts, reqAllocPercent);
                    double allocPerc = actualAlloc / sample.getElevatedCount();

                    if (allocPerc >= reqAllocPercent)
                    {
                        bucketGroup.addSample(sampleId, allocCounts, false);
                        sample.addElevBucketGroup(bucketGroup, allocPerc);
                        wasAllocated = true;
                        ++allocatedSamples;
                        mReassessSamples.add(sampleId);

                        LOGGER.debug(String.format("sample(%d) skipped now added to bg(%d) buckets(grp=%d sam=%d unalloc=%d) count(prop=%s act=%s of %s) allocatedPerc(+%.3f -> %.3f) groupCount(%d)",
                                sampleId, bucketGroup.getId(), bucketGroup.getBucketIds().size(), sample.getElevatedBuckets().size(),
                                sample.getUnallocBuckets().size(), sizeToStr(proposedAllocTotal), sizeToStr(actualAlloc),
                                sizeToStr(sample.getElevatedCount()), sample.lastAllocPercChange(), sample.getAllocPercent(), sample.getElevBucketGroups().size()));
                    }
                }

                if (!wasAllocated)
                {
                    if (maxOtherGroupAllocPerc >= reqAllocPercent)
                    {
                        //  LOGGER.debug(String.format("sample(%d) skipped again with better allocation(%s perc=%.3f)",
                        //      sampleId, sizeToStr(maxOtherGroupAlloc), maxOtherGroupAllocPerc));

                        ++samIndex;
                        continue;
                    }
                    else
                    {
                        LOGGER.debug("sample({}) skipped & not allocated, being removed", sampleId);
                    }
                }

                mSkippedSamples.remove(samIndex);
            }
        }

        // finally add in the newly skipped samples
        for(Integer sampleId : newSkippedSamples)
        {
            if(!mSkippedSamples.contains(sampleId))
                mSkippedSamples.add(sampleId);
        }

        if(initSkippedSamples > mSkippedSamples.size())
        {
            LOGGER.debug("skipped samples: new({}) initial({}) current({}) allocated({}) reskipped({})",
                    newSkippedSamples.size(), initSkippedSamples, mSkippedSamples.size(), allocatedSamples, reskippedSamples);
        }
    }

    private void clearBucketGroups()
    {
        if(mBucketGroups.isEmpty())
            return;

        int initBgCount = mBucketGroups.size();

        int bgIndex = 0;
        while(bgIndex < mBucketGroups.size())
        {
            BucketGroup bucketGroup = mBucketGroups.get(bgIndex);

            if(mFinalBucketGroups.contains(bucketGroup))
            {
                mBucketGroups.remove(bgIndex);
                continue;
            }

            /*
            if(mSkippedBucketGroups.contains(bucketGroup))
            {
                ++bgIndex;
                continue;
            }
            */

            // if any of the founding samples has now been allocated, this group
            // can be purged. Those sample will re-create it during the group-formation routine if warranted
            boolean initSampleAllocated = false;
            for(Integer initSampleId : bucketGroup.getInitialSampleIds())
            {
                if(mReassessSamples.contains(initSampleId))
                {
                    initSampleAllocated = true;
                    break;
                }
            }

            if(initSampleAllocated)
            {
                mBucketGroups.remove(bgIndex);
                continue;
            }

            // remove any samples from this group which were just allocated
            // int sharedCount = 0;
            // double bgSampleCount = bucketGroup.getSampleIds().size();
            for(int samIndex = 0; samIndex < bucketGroup.getSampleIds().size(); ++samIndex)
            {
                int sampleId = bucketGroup.getSampleIds().get(samIndex);

                if(!mReassessSamples.contains(sampleId))
                    continue;

                SampleData sample = mSampleData.get(sampleId);

                // ++sharedCount;

                // remove this sample's allocation to the group
                boolean ok = bucketGroup.removeSampleAllocation(samIndex, sampleId, sample.getElevatedCount());

                if(!ok)
                {
                    LOGGER.debug("bg({}) removal of sample({})", bucketGroup.getId(), sampleId);
                    mHasErrors = true;
                    return;
                }
            }

            // double sharedPerc = sharedCount/bgSampleCount;

            if(bucketGroup.getSampleIds().isEmpty())
            {
                // assume this group needs reassessing since the bulk of its samples have just been allocated
                mBucketGroups.remove(bgIndex);
                continue;
            }

            ++bgIndex;
        }

        if(initBgCount > mBucketGroups.size())
        {
            LOGGER.debug("bucket groups clean-up({} -> {}) from {} reassess samples, skipped groups({})",
                    initBgCount, mBucketGroups.size(), mReassessSamples.size(), mSkippedBucketGroups.size());
        }
    }

    private void applyPredefinedSigs()
    {
        if (mPredefinedSigs == null || mApplyPredefinedSigCount < 1)
            return;

        int bgSigCount = mSpecificCancer.isEmpty() ? mBackgroundGroups.size() : mCancerSamplesMap.size();

        // assume that the first X sigs are background, so start with the first elevated sig
        int sigsApplied = 0;
        for(int sig = bgSigCount; sig < mPredefinedSigs.Cols; ++sig)
        {
            final double[] bucketRatios = mPredefinedSigs.getCol(sig);

            BucketGroup bucketGroup = new BucketGroup(++mNextBucketId);
            bucketGroup.setTag("predefined");
            bucketGroup.setBucketRatios(bucketRatios);

            List<Integer> bucketIds = Lists.newArrayList();

            for(int b = 0; b < mBucketCount; ++b)
            {
                if(bucketRatios[b] == 0)
                    continue;

                bucketIds.add(b);
                bucketGroup.addBucket(b, true);
            }

            LOGGER.debug("created predefined bg({}) with {} buckets", bucketGroup.getId(), bucketIds.size());

            if(!mFinalFitOnly)
            {
                // now test this against every sample in the usual manner
                for (SampleData sample : mSampleData)
                {
                    if (sample.isExcluded())
                        continue;

                    final List<Integer> samBuckets = sample.getElevatedBuckets();

                    double reqAllocPercent = minAllocPercent(sample, false);

                    if (samBuckets.isEmpty() || sample.getUnallocPercent() < reqAllocPercent)
                        continue;

                    double[] allocCounts = sample.getPotentialUnallocCounts(bucketRatios, bucketGroup.getBucketIds());

                    double proposedAllocTotal = sumVector(allocCounts);
                    double proposedAllocPerc = proposedAllocTotal / sample.getElevatedCount();

                    if (proposedAllocPerc < reqAllocPercent || proposedAllocTotal < mMinSampleAllocCount)
                        continue;

                    // apply to the remaining unallocated elevated counts for this sample
                    double actualAlloc = sample.allocateBucketCounts(allocCounts, reqAllocPercent);
                    double allocPerc = actualAlloc / sample.getElevatedCount();

                    if (allocPerc >= reqAllocPercent)
                    {
                        bucketGroup.addSample(sample.Id, allocCounts, false);
                        sample.addElevBucketGroup(bucketGroup, allocPerc);

                        LOGGER.debug(String.format("sample(%d) added to predefined bg(%d) count(prop=%s act=%s of %s) allocatedPerc(+%.3f -> %.3f) groupCount(%d)",
                                sample.Id, bucketGroup.getId(), sizeToStr(proposedAllocTotal), sizeToStr(actualAlloc), sizeToStr(sample.getElevatedCount()),
                                sample.lastAllocPercChange(), sample.getAllocPercent(), sample.getElevBucketGroups().size()));

                        if (!mReassessSamples.contains(sample.Id))
                            mReassessSamples.add(sample.Id);
                    }
                }

                if (bucketGroup.getSampleIds().isEmpty())
                    continue;
            }

            mFinalBucketGroups.add(bucketGroup);

            ++sigsApplied;

            if(sigsApplied >= mApplyPredefinedSigCount)
                break;
        }

        if(!mFinalBucketGroups.isEmpty())
        {
            analyseGroupsVsExtData(mFinalBucketGroups, false);
            logBucketGroups(true);
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

        SigContributionOptimiser sigOptim = new SigContributionOptimiser(mBucketCount, false, SAMPLE_ALLOCATED_PERCENT, true);

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

        List<double[]> ratiosCollection = Lists.newArrayList();
        List<Integer> sigIds = Lists.newArrayList();

        for (SampleData sample : mSampleData)
        {
            if (sample.isExcluded())
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
                double[] allocCounts = sample.getPotentialUnallocCounts(bucketGroup.getBucketRatios(), bucketGroup.getBucketIds());
                double allocTotal = sumVector(allocCounts);

                if (allocTotal / sampleCount < reqAllocPercent)
                    continue;

                if(allocTotal < mMinSampleAllocCount)
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

            ratiosCollection.clear();
            sigIds.clear();
            double[] potentialGroupContribs = new double[groupCount];

            final double[] sampleCounts = sample.getElevatedBucketCounts();
            final double[] countsNoise = sample.getCountRanges();
            int backgroundGroupIndex = -1;

            for (index = 0; index < groupCount; ++index)
            {
                int bgIndex = bgIndexList.get(index);
                final BucketGroup bucketGroup = mFinalBucketGroups.get(bgIndex);

                ratiosCollection.add(bucketGroup.getBucketRatios());
                sigIds.add(bucketGroup.getId());

                potentialGroupContribs[index] = potentialAllocTotals.get(index);

                if(bucketGroup.equals(sample.getBackgroundGroup()))
                    backgroundGroupIndex = index;
            }

            if(mSampleWatchList.contains(sample.Id))
                sigOptim.setLogVerbose(true);
            else
                sigOptim.setLogVerbose(false);

            sigOptim.initialise(sample.Id, sampleCounts, countsNoise, ratiosCollection, potentialGroupContribs, MIN_GROUP_ALLOC_PERCENT_LOWER, mMinSampleAllocCount);
            sigOptim.setSigIds(sigIds);

            // each sample's background sig will remain in the list even if it drops below the required threshold
            sigOptim.setRequiredSig(backgroundGroupIndex);
            boolean validCalc = sigOptim.fitToSample();

            if(!validCalc)
            {
                LOGGER.error("sample({}) refit of {} sigs failed", sample.Id, groupCount);
                mHasErrors = true;
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
            double maxGroupPerc = 0;

            for(index = 0; index < sortedAllocIndices.size(); ++index)
            {
                int itemIndex = sortedAllocIndices.get(index);
                Integer bgIndex = bgIndexList.get(itemIndex);
                final BucketGroup bucketGroup = mFinalBucketGroups.get(bgIndex);
                double potentialAlloc = potentialGroupContribs[itemIndex];
                double fitAlloc = newGroupContribs[itemIndex];

                if(fitAlloc/sampleCount < reqAllocPercent || fitAlloc < mMinSampleAllocCount)
                {
                    if(useFittedAllocs && fitAlloc > 0)
                    {
                        LOGGER.warn(String.format("sample(%d) missed fit contrib bg(%d) count(pot=%s fit=%s perc=%.3f of %s)",
                                sample.Id, bucketGroup.getId(), sizeToStr(potentialAlloc), sizeToStr(fitAlloc), fitAlloc / sampleCount, sizeToStr(sampleCount)));
                    }
                    break;
                }

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

                double actualAlloc = sample.allocateBucketCounts(allocCounts, reqAllocPercent);
                double allocPerc = actualAlloc / sample.getElevatedCount();

                if(allocPerc >= reqAllocPercent)
                {
                    maxGroupPerc = max(allocPerc, maxGroupPerc);
                    bucketGroup.addSample(sample.Id, allocCounts, false);
                    sample.addElevBucketGroup(bucketGroup, allocPerc);

                    LOGGER.debug(String.format("sample(%d) added to bg(%d) count(pot=%s fit=%s act=%s of %s) allocatedPerc(+%.3f -> %.3f) noise(%s %.3f/%.3f) groupCount(%d)",
                            sample.Id, bucketGroup.getId(), sizeToStr(potentialAlloc), sizeToStr(fitAlloc), sizeToStr(actualAlloc), sizeToStr(sampleCount),
                            sample.lastAllocPercChange(), sample.getAllocPercent(), sizeToStr(sample.getAllocNoise()), sample.getNoisePerc(), sample.getNoiseOfTotal(),
                            sample.getElevBucketGroups().size()));
                }
            }

            String allocPercChange = "unch";

            if(sample.getAllocPercent() > prevTotalAllocPerc + 0.01)
                allocPercChange = "better";
            else if(sample.getAllocPercent() < prevTotalAllocPerc - 0.01)
                allocPercChange = "worse";

            LOGGER.debug(String.format("sample(%d) final fit: method(%s) groups(%d prev=%d max=%d) %s allocation(prev=%.3f new=%.3f, act=%s of %s) proposed(%s) fit(%s)",
                    sample.Id, useFittedAllocs ? "fit" : "best", sample.getElevBucketGroups().size(), prevGroupCount, groupCount,
                    allocPercChange, prevTotalAllocPerc, sample.getAllocPercent(), sizeToStr(sample.getAllocatedCount()),
                    sizeToStr(sample.getElevatedCount()), sizeToStr(maxAllocTotal), sizeToStr(fitAllocTotal)));

            // log any missed allocations
            /*
            if(sample.getAllocPercent() < 0.9)
            {
                for (final BucketGroup bucketGroup : mSkippedBucketGroups)
                {
                    if (!bucketGroup.hasSample(sample.Id))
                        continue;

                    double skippedAlloc = bucketGroup.getSampleCount(sample.Id);
                    double skippedAllocPerc = skippedAlloc/sample.getElevatedCount();

                    if (skippedAllocPerc > MIN_GROUP_ALLOC_PERCENT && skippedAlloc > sample.getUnallocPercent())
                    {
                        LOGGER.debug(String.format("sample(%d) missed skipped bg(%d) alloc(%s perc=%.3f) vs actualMax(%.3f)",
                                sample.Id, bucketGroup.getId(), sizeToStr(skippedAlloc), skippedAllocPerc, maxGroupPerc));
                    }
                }
            }
            */
        }

        // LOGGER.debug(String.format("small alloc total(%s)", sizeToStr(smallAllocsTotal)));
    }

    private double minAllocPercent(final SampleData sample, boolean scaleToLower)
    {
        if(!scaleToLower)
            return MIN_GROUP_ALLOC_PERCENT;

        // require say 10% of the remaining unallocated count, but with an upper minimum limit
        return max(sample.getUnallocPercent() * MIN_GROUP_ALLOC_PERCENT, MIN_GROUP_ALLOC_PERCENT_LOWER);
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

    private void analyseGroupsVsExtData(List<BucketGroup> bgList, boolean verbose)
    {
        if(mExtSampleData == null)
            return;

        LOGGER.debug("analysing {} bucket groups", bgList.size());

        // check counts for each external data category against the samples in each bucket group

        // want to know if:
        // a) the samples have 1 or more identical categories and
        // b) the proportion of these categories out of the cohort (ie missing other samples, linked to other bucket groups)

        for (BucketGroup bucketGroup : bgList)
        {
            final List<Integer> sampleIds = bucketGroup.getSampleIds();
            int sampleCount = sampleIds.size();

            // clear any previous
            bucketGroup.setCancerType("");
            bucketGroup.setEffects("");

            String groupEffectsStr = "";
            String groupCancerTypeStr = "";
            String catEffectsStr = "";
            String catCancerTypeStr = "";

            HashMap<String,Integer> categoryCounts = populateSampleCategoryMap(sampleIds);

            for(Map.Entry<String,Integer> entrySet : categoryCounts.entrySet())
            {
                final String catName = entrySet.getKey();
                int catSampleCount = entrySet.getValue();
                double samplesPerc = catSampleCount/(double)sampleCount; // percentage of samples within the group

                int catAllCount = mExtCategoriesMap.get(catName);
                double catPerc = catSampleCount / (double) catAllCount; // percentage of all instances of this

                // ignore small & rare effects unless they're the majority of this group
                boolean ignoreCatPercent = catAllCount < 0.5 * sampleCount && catAllCount < 30;
                boolean hasDominantSamplePerc = samplesPerc >= DOMINANT_CATEGORY_PERCENT;
                boolean hasProminentSamplePerc = samplesPerc >= DOMINANT_CATEGORY_PERCENT * 0.5;
                boolean hasDominantCatPerc = catPerc >= DOMINANT_CATEGORY_PERCENT && !ignoreCatPercent;
                boolean hasProminentCatPerc = (catPerc >= DOMINANT_CATEGORY_PERCENT * 0.5) && !ignoreCatPercent;

                if(!hasProminentSamplePerc && !hasProminentCatPerc)
                    continue;

                if(verbose)
                {
                    LOGGER.debug(String.format("bg(%d) category(%s) count(%d of %d) perc(group=%.3f category=%.3f)",
                            bucketGroup.getId(), catName, catSampleCount, sampleCount, samplesPerc, catPerc));
                }

                if(extractCategoryName(catName).equals(CATEGORY_CANCER_TYPE))
                {
                    String cancerType = extractCategoryValue(catName);

                    if(hasDominantSamplePerc)
                    {
                        groupCancerTypeStr = String.format("%s=%.2f", cancerType, samplesPerc);
                    }
                    else if(hasProminentSamplePerc)
                    {
                        if(!groupCancerTypeStr.isEmpty())
                            groupCancerTypeStr += ";";

                        groupCancerTypeStr += String.format("%s=%.2f", cancerType, samplesPerc);
                    }

                    if(hasDominantCatPerc)
                    {
                        if(!catCancerTypeStr.isEmpty())
                            catCancerTypeStr += ";";

                        catCancerTypeStr += String.format("%s=%.2f", cancerType, catPerc);
                    }
                }
                else
                {
                    String effect = catName.length() > 10 ? catName.substring(0, 10) : catName;

                    if(hasDominantSamplePerc)
                    {
                        if (!groupEffectsStr.isEmpty())
                            groupEffectsStr += ";";

                        groupEffectsStr += String.format("%s=%.2f", effect, samplesPerc);
                    }

                    if(hasDominantCatPerc)
                    {
                        if(!catEffectsStr.isEmpty())
                            catEffectsStr += ";";

                        catEffectsStr += String.format("%s=%.2f", effect, catPerc);
                    }
                }
            }

            String cancerTypeStr = groupCancerTypeStr;

            if(!catCancerTypeStr.isEmpty())
            {
                if(!cancerTypeStr.isEmpty())
                    cancerTypeStr += ";";

                cancerTypeStr += "cat:" + catCancerTypeStr;
            }

            bucketGroup.setCancerType(cancerTypeStr);

            String effectsStr = groupEffectsStr;

            if(!catEffectsStr.isEmpty())
            {
                if(!effectsStr.isEmpty())
                    effectsStr += ";";

                effectsStr += "cat:" + catEffectsStr;
            }

            bucketGroup.setEffects(effectsStr);
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

                LOGGER.debug(String.format("rank %d: bg(%d) %scancer(%s) samples(%d) variants(avg=%s total=%s perc=%.3f) buckets(%d: %s) effects(%s)",
                        i, bucketGroup.getId(), linkData, bucketGroup.getCancerType(), bucketGroup.getSampleIds().size(),
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

    private void writeInterimBucketGroups()
    {
        try
        {
            if(mBucketGroupFileWriter == null)
            {
                mBucketGroupFileWriter = getNewFile(mOutputDir, mOutputFileId + "_ba_interim_groups.csv");

                mBucketGroupFileWriter.write("RunId,Rank,BgId,Type,CancerType,Effects,SampleCount,BucketCount,MutLoad,PotAlloc,PotAllocAdj");

                for(int i = 0; i < mBucketCount; ++i)
                {
                    mBucketGroupFileWriter.write(String.format(",%d", i));
                }

                mBucketGroupFileWriter.newLine();
            }

            BufferedWriter writer = mBucketGroupFileWriter;

            for (int bgIndex = 0; bgIndex < mTopAllocBucketGroups.size(); ++bgIndex)
            {
                BucketGroup bucketGroup = mTopAllocBucketGroups.get(bgIndex);

                String groupType = bucketGroup.getTag().equalsIgnoreCase("background") ? "BG" : "Elev";

                writer.write(String.format("%d,%d,%d,%s,%s,%s,%d,%d,%.0f,%.0f,%.0f",
                        mRunId, bgIndex, bucketGroup.getId(), groupType, bucketGroup.getCancerType(), bucketGroup.getEffects(),
                        bucketGroup.getSampleIds().size(), bucketGroup.getBucketIds().size(),
                        sumVector(bucketGroup.getBucketCounts()), bucketGroup.getPotentialAllocation(), bucketGroup.getPotentialAdjAllocation()));

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
            LOGGER.error("failed to write output file: interim bucket groups");
        }
    }

    private void writeFinalBucketGroups()
    {
        try
        {
            BufferedWriter writer = getNewFile(mOutputDir, mOutputFileId + "_ba_group_data.csv");

            writer.write("Rank,BgId,Type,CancerType,Effects,SampleCount,BucketCount,MutLoad,PotentialAlloc,Purity");

            for(int i = 0; i < mBucketCount; ++i)
            {
                writer.write(String.format(",%d", i));
            }

            writer.newLine();

            for (int bgIndex = 0; bgIndex < mFinalBucketGroups.size(); ++bgIndex)
            {
                BucketGroup bucketGroup = mFinalBucketGroups.get(bgIndex);

                String groupType = bucketGroup.getTag().equalsIgnoreCase("background") ? "BG" : "Elev";

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

            writer.close();

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
        double noisAllocated = 0;

        for(final SampleData sample : mSampleData)
        {
            mAllocatedCount += sample.getAllocatedCount();

            if(sample.getAllocPercent() >= SAMPLE_ALLOCATED_PERCENT)
                ++fullyAllocated;

            noisAllocated += sample.getAllocNoise();
        }

        LOGGER.debug(String.format("overall: samples(%d alloc=%d) groups(%d) counts: total(%s) background(%s perc=%.3f) elevated(%s perc=%.3f) allocated(%s perc=%.3f) allocNoise(%s perc=%.3f)",
                mSampleCount, fullyAllocated, mFinalBucketGroups.size(), sizeToStr(mTotalCount), sizeToStr(mBackgroundCount), mBackgroundCount/mTotalCount,
                sizeToStr(mElevatedCount), mElevatedCount/mTotalCount, sizeToStr(mAllocatedCount), mAllocatedCount/mElevatedCount,
                sizeToStr(noisAllocated), noisAllocated/mElevatedCount));
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

        LOGGER.debug(String.format("sample summary: total(%d) elev(%d) alloc(%d, %.3f of %s) partial(%d) unalloc(%d)",
                mSampleCount, mAllSampleBucketGroups.size(), fullyAllocCount, percAllocated, sizeToStr(mElevatedCount),
                partiallyAllocCount, noMatchCount));

        logOverallStats();
    }

    private void logWorstAllocatedSamples()
    {
        List<int[]> worstAllocatedSamples = Lists.newArrayList();
        double[] unallocatedByBucket = new double[mBucketCount];

        int SAMPLE_ID = 0;
        int UNALLOC_AMT = 1;

        int nonExcludedSampleCount = 0;

        for(final SampleData sample : mSampleData)
        {
            if(sample.isExcluded())
                continue;

            ++nonExcludedSampleCount;
            sumVectors(sample.getUnallocBucketCounts(), unallocatedByBucket);

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

        LOGGER.debug(String.format("worst %d (perc=%.3f) samples have unallocted(%s perc=%.3f of %s)",
                sampleCount, sampleCount/(double)nonExcludedSampleCount, sizeToStr(totalUnallocated), totalUnallocated/mElevatedCount, sizeToStr(mElevatedCount)));

        // report worst allocated buckets in a similar way
        List<Integer> worstBucketIndices = getSortedVectorIndices(unallocatedByBucket, false);
        for(int worstIndex = 0; worstIndex < 10; ++worstIndex)
        {
            int worstBucket = worstBucketIndices.get(worstIndex);
            double bucketTotal = unallocatedByBucket[worstBucket];
            double totalBucketCount = sumVector(mElevatedCounts.getRow(worstBucket));

            LOGGER.debug(String.format("%d: worst bucket(%d) unallocated(%s of %s, perc=%.3f) percOfTotal(%.4f)",
                    worstIndex, worstBucket, sizeToStr(bucketTotal), sizeToStr(totalBucketCount), bucketTotal/totalBucketCount, bucketTotal/mElevatedCount));
        }
    }

    private void logSimilarSampleContribs()
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

        // exclude samples with onl
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

                List<Integer> sam1BgList = sample1.getElevBucketGroups().stream().map(BucketGroup::getId).collect(Collectors.toList());
                List<Integer> sam2BgList = sample2.getElevBucketGroups().stream().map(BucketGroup::getId).collect(Collectors.toList());

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

                if(i < 20)
                {
                    LOGGER.debug(String.format("samples(%d and %d) similar sig contribs css(%.4f) groups(s1=%d s2=%d match=%d) ids(%s)",
                            samId1, samId2, result[CSSR_VAL], sam1BgList.size(), sam2BgList.size(), bgOverlaps.size(), bgIdsStr));
                }
            }

            for(Map.Entry<String,Integer> entry : groupComboMap.entrySet())
            {
                if(entry.getValue() < 5)
                    continue;

                LOGGER.debug("bgGroups({}) repeated {} times", entry.getKey(), entry.getValue());
            }
        }
    }

    private void compareSignatures()
    {
        double sigCompareCss = 0.95;

        if(mProposedSigs != null)
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
                    int sigId1 = (int)result[CSSR_I1];
                    int sigId2 = (int)result[CSSR_I2];
                    final BucketGroup bg1 = mFinalBucketGroups.get(mSigToBgMapping.get(sigId1));
                    final BucketGroup bg2 = mFinalBucketGroups.get(mSigToBgMapping.get(sigId2));

                    // ignore reporting similarities in the BG groups
                    if(mBackgroundGroups.contains(bg1) && mBackgroundGroups.contains(bg2))
                        continue;

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

    private void createSignatures()
    {
        if(mMaxProposedSigs == 0)
            return;

        int proposedSigCount = min(mMaxProposedSigs, mFinalBucketGroups.size());

        LOGGER.debug("creating {} signatures", proposedSigCount);

        mProposedSigs = new NmfMatrix(mBucketCount, proposedSigCount);

//        NmfMatrix sigTest = new NmfMatrix(mBucketCount, 1);
//
//        double purityThreshold = 0.7;
//        int scoreThreshold = 50; // eg 100% pure, 2 samples 10 buckets, LF=2.5
//        int minSamples = 2;

        int sigId = 0;

        List<String> sigNames = Lists.newArrayList();

        for(int bgIndex = 0; bgIndex < mFinalBucketGroups.size(); ++bgIndex)
        {
            final BucketGroup bucketGroup = mFinalBucketGroups.get(bgIndex);

            String sigName = "";

            if(mBackgroundGroups.contains(bucketGroup))
            {
                sigName = "BG_";

                if(bucketGroup.getCancerType().length() > 5)
                    sigName += bucketGroup.getCancerType().substring(0,5);
                else
                    sigName += bucketGroup.getCancerType();

                sigName = sigName.replaceAll("/", "");
            }
            else
            {
                // take any dominant effect
                sigName = String.format("ELV_%d", bucketGroup.getId());

                String details = "";
                if(!bucketGroup.getCancerType().isEmpty())
                {
                    details += "_" + bucketGroup.getCancerType().substring(0,3);
                }

                if(!bucketGroup.getEffects().isEmpty())
                {
                    if(bucketGroup.getEffects().length() > 5)
                        details += "_" + bucketGroup.getEffects().substring(0,5);
                    else
                        details += "_" + bucketGroup.getEffects();

                }

                sigName += details.replaceAll("[0123456789.=]", "");
            }

            sigNames.add(sigName);


            final double[] bucketRatios = bucketGroup.getBucketRatios();

            /*
            if(bucketGroup.getPurity() < purityThreshold || bucketGroup.calcScore() < scoreThreshold
            || bucketGroup.getSampleIds().size() < minSamples)
                continue;


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

            */

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

        writeSignatures(mProposedSigs, "_ba_sigs.csv", sigNames);
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
        if (mSampleData.isEmpty())
            return;

        mCancerSamplesMap = new HashMap();

        List<String> cancerTypes = Lists.newArrayList();

        for (final SampleData sample : mSampleData)
        {
            final String cancerType = sample.getCancerType();
            List<Integer> samplesList = mCancerSamplesMap.get(cancerType);

            if (samplesList == null)
            {
                samplesList = Lists.newArrayList();
                mCancerSamplesMap.put(cancerType, samplesList);
                cancerTypes.add(cancerType);
            }

            samplesList.add(sample.Id);
        }

        // put all low sample count by cancer type into the 'Other' cancer type mapping group
        String otherType = "Other";
        int minSamples = 20;

        List<Integer> otherTypeList = mCancerSamplesMap.get(otherType);
        if (otherTypeList == null)
        {
            otherTypeList = Lists.newArrayList();
        }

        for (final String cancerType : cancerTypes)
        {
            if (cancerType.equals(otherType))
                continue;

            List<Integer> samplesList = mCancerSamplesMap.get(cancerType);

            if (samplesList.size() < minSamples)
            {
                otherTypeList.addAll(samplesList);
                mCancerSamplesMap.remove(cancerType);
            }
        }

        if (!otherTypeList.isEmpty())
        {
            mCancerSamplesMap.put(otherType, otherTypeList);
        }

        // check all samples were allocated
        int sampleCount = 0;
        for (Map.Entry<String, List<Integer>> entry : mCancerSamplesMap.entrySet())
        {
            final List<Integer> sampleIds = entry.getValue();
            sampleCount += sampleIds.size();
        }

        if(sampleCount != mSampleCount)
        {
            LOGGER.error("sample count mismatch({} vs {})", sampleCount, mSampleCount);
            mHasErrors = true;
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
                    continue;

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
                // write some meaningful names for the sigs
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

    public void writeSampleContributions()
    {
        int sigCount = mFinalBucketGroups.size();

        NmfMatrix contribMatrix = new NmfMatrix(sigCount, mSampleCount);
        double[][] contribData = contribMatrix.getData();

        for(int bgIndex = 0; bgIndex < mFinalBucketGroups.size(); ++bgIndex)
        {
            final BucketGroup bucketGroup = mFinalBucketGroups.get(bgIndex);
            final List<Integer> sampleIds = bucketGroup.getSampleIds();
            final List<Double> sampleContribs = bucketGroup.getSampleCountTotals();

            for(int samIndex = 0; samIndex < sampleIds.size(); ++samIndex)
            {
                Integer sampleId = sampleIds.get(samIndex);
                double sampleContrib = sampleContribs.get(samIndex);

                contribData[bgIndex][sampleId] = sampleContrib;
            }
        }

        contribMatrix.cacheTranspose();

        // scale so the contribs do not exceed the actual sample count totals
        /*
        for(int s = 0; s < mSampleCount; ++s)
        {
            double adjRatio = 0;
            double sampleTotal = mSampleTotals[s];
            double contribTotal = sumVector(contribMatrix.getCol(s));

            if(contribTotal > sampleTotal)
            {
                adjRatio = sampleTotal / contribTotal;

                for(int sig = 0; sig < sigCount; ++sig)
                {
                    contribData[sig][s] *= adjRatio;
                }
            }
        }
        */

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

            writeMatrixData(writer, contribMatrix, false);

            writer.close();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing sample contrib file");
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

                writer.newLine();
            }

            writer.close();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile");
        }
    }


}
