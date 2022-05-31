package com.hartwig.hmftools.cup.somatics;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;
import static com.hartwig.hmftools.common.stats.CosineSimilarity.calcCosineSim;
import static com.hartwig.hmftools.common.stats.Percentiles.getPercentile;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.VariantType.SNP;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaRefFiles.purpleSomaticVcfFile;
import static com.hartwig.hmftools.cup.common.CategoryType.SAMPLE_TRAIT;
import static com.hartwig.hmftools.cup.common.CategoryType.SNV;
import static com.hartwig.hmftools.cup.common.ClassifierType.GENOMIC_POSITION_COHORT;
import static com.hartwig.hmftools.cup.common.ClassifierType.SNV_96_PAIRWISE;
import static com.hartwig.hmftools.cup.common.CupCalcs.adjustRefCounts;
import static com.hartwig.hmftools.cup.common.CupCalcs.calcPercentilePrevalence;
import static com.hartwig.hmftools.cup.common.CupCalcs.convertToPercentages;
import static com.hartwig.hmftools.cup.common.CupConstants.CSS_SIMILARITY_CUTOFF;
import static com.hartwig.hmftools.cup.common.CupConstants.CSS_SIMILARITY_MAX_MATCHES;
import static com.hartwig.hmftools.cup.common.CupConstants.DATA_TYPE_SNV_COUNT;
import static com.hartwig.hmftools.cup.common.CupConstants.GEN_POS_BUCKET_SIZE;
import static com.hartwig.hmftools.cup.common.CupConstants.GEN_POS_CSS_EXPONENT_CUTOFF;
import static com.hartwig.hmftools.cup.common.CupConstants.GEN_POS_CSS_EXPONENT_TAIL;
import static com.hartwig.hmftools.cup.common.CupConstants.GEN_POS_MAX_SAMPLE_COUNT;
import static com.hartwig.hmftools.cup.common.CupConstants.SNV_96_CSS_DIFF_EXPONENT;
import static com.hartwig.hmftools.cup.common.CupConstants.SNV_96_CSS_THRESHOLD;
import static com.hartwig.hmftools.cup.common.CupConstants.GEN_POS_CSS_THRESHOLD;
import static com.hartwig.hmftools.cup.common.CupConstants.GEN_POS_CSS_EXPONENT;
import static com.hartwig.hmftools.cup.common.CupConstants.UNDEFINED_PERC_MAX_MULTIPLE;
import static com.hartwig.hmftools.cup.common.ResultType.CLASSIFIER;
import static com.hartwig.hmftools.cup.common.ResultType.LIKELIHOOD;
import static com.hartwig.hmftools.cup.common.ResultType.PERCENTILE;
import static com.hartwig.hmftools.cup.common.SampleData.isKnownCancerType;
import static com.hartwig.hmftools.cup.common.SampleResult.checkIsValidCancerType;
import static com.hartwig.hmftools.cup.common.SampleSimilarity.recordCssSimilarity;
import static com.hartwig.hmftools.cup.somatics.GenomicPositions.convertSomaticVariantsToPosFrequencies;
import static com.hartwig.hmftools.cup.somatics.SomaticDataLoader.loadRefSampleCounts;
import static com.hartwig.hmftools.cup.somatics.SomaticDataLoader.loadRefSignaturePercentileData;
import static com.hartwig.hmftools.cup.somatics.SomaticDataLoader.loadSampleCountsFromFile;
import static com.hartwig.hmftools.cup.somatics.SomaticDataLoader.loadSampleMatrixData;
import static com.hartwig.hmftools.cup.somatics.SomaticDataLoader.loadSomaticVariants;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.INCLUDE_AID_APOBEC_SIG;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.INCLUDE_AID_APOBEC_SIG_DESC;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.INCLUDE_AID_APOBEC;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.INCLUDE_AID_APOBEC_DESC;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.applyMaxCssAdjustment;
import static com.hartwig.hmftools.cup.somatics.TrinucleotideCounts.convertSomaticVariantsToSnvCounts;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sigs.PositionFrequencies;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.cup.CuppaConfig;
import com.hartwig.hmftools.cup.common.CategoryType;
import com.hartwig.hmftools.cup.common.CuppaClassifier;
import com.hartwig.hmftools.cup.common.NoiseRefCache;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;
import com.hartwig.hmftools.cup.common.SampleSimilarity;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class SomaticClassifier implements CuppaClassifier
{
    private final CuppaConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private Matrix mRefSampleSnv96Counts;
    private final Map<String,Integer> mRefSampleSnv96CountsIndex;
    private final Map<String,double[]> mRefCancerSnvCountPercentiles;

    private Matrix mRefCancerGenPosCounts;
    private final List<String> mRefGenPosCancerTypes;

    private Matrix mRefSampleGenPosCounts;
    private final Map<String,Integer> mRefSampleGenPosCountsIndex;

    private Matrix mSampleSnv96Counts;
    private final Map<String,Integer> mSampleSnv96CountsIndex; // index of a sampleId into the sample SNV counts matrix
    private final Map<String,Integer> mSampleSnvTotals;

    private Matrix mSampleGenPosCounts;
    private final Map<String,Integer> mSampleGenPosCountsIndex;

    private final SigContributions mSigContributions;
    private final PositionFrequencies mPosFrequencies;

    private boolean mIsValid;
    private BufferedWriter mGenPosCohortCssWriter;

    private final boolean mIncludeAidApobecGenPos;

    private final boolean mRunPairwiseGenPos;
    private final double mMaxCssAdjustFactorSnv;
    private final double mMaxCssAdjustFactorGenPos;
    private final double mCssExponentSnv;
    private final double mGenPosCssExponent;
    private final double mGenPosCssExponentTail;
    private final double mGenPosCssExponentCutoff;
    private final boolean mWriteGenPosSims;
    private final boolean mWriteSnvSims;

    // config
    public static final String MAX_CSS_ADJUST_FACTOR_SNV = "css_max_factor_snv";
    public static final String MAX_CSS_ADJUST_FACTOR_GEN_POS = "css_max_factor_gen_pos";

    public static final String GEN_POS_PAIRWISE = "gen_pos_pairwise";
    public static final String CSS_EXPONENT_SNV = "snv_css_exponent";
    public static final String CSS_EXPONENT_GEN_POS = "gen_pos_css_exponent";
    public static final String CSS_EXPONENT_GEN_POS_TAIL = "gen_pos_css_exponent_tail";
    public static final String CSS_EXPONENT_GEN_POS_CUTOFF = "gen_pos_css_exponent_cutoff";

    public static final String SNV_POS_FREQ_POS_SIZE = "gen_pos_bucket_size";

    public static final String WRITE_GEN_POS_CSS = "write_gen_pos_css";
    public static final String WRITE_GEN_POS_SIMILARITIES = "write_gen_pos_sims";
    public static final String WRITE_SNV_SIMILARITIES = "write_snv_sims";

    public SomaticClassifier(final CuppaConfig config, final SampleDataCache sampleDataCache, final CommandLine cmd)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mSampleSnv96Counts = null;
        mRefCancerGenPosCounts = null;
        mRefSampleGenPosCounts = null;
        mRefSampleSnv96CountsIndex = Maps.newHashMap();
        mSampleSnv96CountsIndex = Maps.newHashMap();
        mSampleGenPosCountsIndex = Maps.newHashMap();
        mSampleSnvTotals = Maps.newHashMap();

        mRefSampleSnv96Counts = null;
        mRefCancerSnvCountPercentiles = Maps.newHashMap();
        mRefGenPosCancerTypes = Lists.newArrayList();
        mRefSampleGenPosCountsIndex = Maps.newHashMap();

        if(cmd != null)
        {
            mIncludeAidApobecGenPos = cmd.hasOption(INCLUDE_AID_APOBEC);
            mRunPairwiseGenPos = cmd.hasOption(GEN_POS_PAIRWISE);

            mCssExponentSnv = Double.parseDouble(cmd.getOptionValue(CSS_EXPONENT_SNV, String.valueOf(SNV_96_CSS_DIFF_EXPONENT)));
            mGenPosCssExponent = Double.parseDouble(cmd.getOptionValue(CSS_EXPONENT_GEN_POS, String.valueOf(GEN_POS_CSS_EXPONENT)));
            mGenPosCssExponentTail = Double.parseDouble(cmd.getOptionValue(CSS_EXPONENT_GEN_POS_TAIL, String.valueOf(GEN_POS_CSS_EXPONENT_TAIL)));
            mGenPosCssExponentCutoff = Double.parseDouble(cmd.getOptionValue(CSS_EXPONENT_GEN_POS_CUTOFF, String.valueOf(GEN_POS_CSS_EXPONENT_CUTOFF)));

            mMaxCssAdjustFactorSnv = Double.parseDouble(cmd.getOptionValue(MAX_CSS_ADJUST_FACTOR_SNV, "0"));
            mMaxCssAdjustFactorGenPos = Double.parseDouble(cmd.getOptionValue(MAX_CSS_ADJUST_FACTOR_GEN_POS, "0"));
            mWriteSnvSims = cmd.hasOption(WRITE_SNV_SIMILARITIES);
            mWriteGenPosSims = cmd.hasOption(WRITE_GEN_POS_SIMILARITIES);
        }
        else
        {
            mIncludeAidApobecGenPos = false;
            mRunPairwiseGenPos = false;
            mCssExponentSnv = SNV_96_CSS_DIFF_EXPONENT;
            mGenPosCssExponent = GEN_POS_CSS_EXPONENT;
            mGenPosCssExponentTail = 0;
            mGenPosCssExponentCutoff = 0;
            mMaxCssAdjustFactorSnv = 0;
            mMaxCssAdjustFactorGenPos = 0;
            mWriteSnvSims = false;
            mWriteGenPosSims = false;
        }

        mIsValid = true;
        mGenPosCohortCssWriter = null;

        mSigContributions = new SigContributions(mConfig, mSampleDataCache);

        int posFreqBucketSize = cmd != null && cmd.hasOption(SNV_POS_FREQ_POS_SIZE) ?
                Integer.parseInt(cmd.getOptionValue(SNV_POS_FREQ_POS_SIZE)) : GEN_POS_BUCKET_SIZE;

        mPosFrequencies = new PositionFrequencies(posFreqBucketSize, GEN_POS_MAX_SAMPLE_COUNT);

        if(mConfig.RefSnvCountsFile.isEmpty() && mConfig.RefSigContributionFile.isEmpty() && mConfig.RefSnvCancerPosFreqFile.isEmpty())
            return;

        loadRefSignaturePercentileData(
                mConfig.RefSigContributionFile, mSigContributions.getRefCancerSigContribPercentiles(), mRefCancerSnvCountPercentiles);

        mRefSampleSnv96Counts = loadSampleCountsFromFile(mConfig.RefSnvCountsFile, mRefSampleSnv96CountsIndex);

        mRefCancerGenPosCounts = loadRefSampleCounts(mConfig.RefSnvCancerPosFreqFile, mRefGenPosCancerTypes, Lists.newArrayList());
        mRefSampleGenPosCounts = loadSampleMatrixData(mConfig.RefSnvSamplePosFreqFile, mRefSampleGenPosCountsIndex);

        if(mRefSampleSnv96Counts == null || mRefCancerGenPosCounts == null)
        {
            CUP_LOGGER.error("invalid somatic matrix data: SNV-96({}) GenPos({})",
                    mRefSampleSnv96Counts != null, mRefCancerGenPosCounts != null);
            mIsValid = false;
        }

        mIsValid &= loadSampleCounts();
        mIsValid &= mSigContributions.loadSigContributions(mSampleSnv96Counts);

        if(mIsValid)
        {
            // record SNV totals, prior to any noise adjustments
            for(Map.Entry<String,Integer> entry : mRefSampleSnv96CountsIndex.entrySet())
            {
                String sampleId = entry.getKey();
                final double[] sampleCounts = mRefSampleSnv96Counts.getRow(entry.getValue());
                int snvTotal = (int) sumVector(sampleCounts);
                mSampleSnvTotals.put(sampleId, snvTotal);
            }

            if(!mConfig.SampleSnvCountsFile.equals(mConfig.RefSnvCountsFile))
            {
                for(Map.Entry<String, Integer> entry : mSampleSnv96CountsIndex.entrySet())
                {
                    String sampleId = entry.getKey();

                    if(mSampleSnvTotals.containsKey(sampleId))
                        continue;

                    final double[] sampleCounts = mSampleSnv96Counts.getRow(entry.getValue());
                    int snvTotal = (int) sumVector(sampleCounts);
                    mSampleSnvTotals.put(sampleId, snvTotal);
                }
            }

            // apply any specified noise
            if(mConfig.NoiseAdjustments.makeNoiseAdjustment(SNV_96_PAIRWISE))
            {
                final double[] noiseAdjustments = mConfig.NoiseAdjustments.getNoiseData(SNV_96_PAIRWISE);
                int noiseAllocation = mConfig.NoiseAdjustments.getNoiseAllocation(SNV_96_PAIRWISE);

                CUP_LOGGER.info("applying noise({}) to SNV-96 counts", noiseAllocation);

                NoiseRefCache.applyNoise(mRefSampleSnv96Counts, noiseAdjustments, noiseAllocation);

                if(mSampleSnv96Counts != mRefSampleSnv96Counts)
                    NoiseRefCache.applyNoise(mSampleSnv96Counts, noiseAdjustments, noiseAllocation);
            }

            if(mRunPairwiseGenPos && mConfig.NoiseAdjustments.makeNoiseAdjustment(GENOMIC_POSITION_COHORT))
            {
                final double[] noiseAdjustments = mConfig.NoiseAdjustments.getNoiseData(GENOMIC_POSITION_COHORT);
                int noiseAllocation = mConfig.NoiseAdjustments.getNoiseAllocation(GENOMIC_POSITION_COHORT);

                CUP_LOGGER.debug("applying noise({}) to genomic position sample counts", noiseAllocation);

                NoiseRefCache.applyNoise(mRefSampleGenPosCounts, noiseAdjustments, noiseAllocation);
            }
        }

        if(cmd.hasOption(WRITE_GEN_POS_CSS) && !mRunPairwiseGenPos)
        {
            initialiseOutputFiles();
        }
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(INCLUDE_AID_APOBEC, false, INCLUDE_AID_APOBEC_DESC);
        options.addOption(GEN_POS_PAIRWISE, false, "Run genomic position as a pairwise classifier");
        options.addOption(INCLUDE_AID_APOBEC_SIG, false, INCLUDE_AID_APOBEC_SIG_DESC);
        options.addOption(MAX_CSS_ADJUST_FACTOR_SNV, true, "Max CSS adustment factor for SNV 96");
        options.addOption(MAX_CSS_ADJUST_FACTOR_GEN_POS, true, "Max CSS adustment factor for genomic pos frequency");
        options.addOption(CSS_EXPONENT_SNV, true, "SNV 96 CSS exponent for tail");
        options.addOption(CSS_EXPONENT_GEN_POS, true, "Genomic position CSS exponent");
        options.addOption(CSS_EXPONENT_GEN_POS_TAIL, true, "Genomic position CSS exponent for tail");
        options.addOption(CSS_EXPONENT_GEN_POS_CUTOFF, true, "Genomic position CSS exponent cutoff");
        options.addOption(SNV_POS_FREQ_POS_SIZE, true, "Genomic position bucket size (default: 20000)");
        options.addOption(WRITE_SNV_SIMILARITIES, false, "Write SNV-96 CSS to file");
        options.addOption(WRITE_GEN_POS_SIMILARITIES, false, "Write genomic position CSS to file");
        options.addOption(WRITE_GEN_POS_CSS, false, "Write gen-pos CSS to file");
    }

    public CategoryType categoryType() { return SNV; }
    public boolean isValid() { return mIsValid; }

    public void close()
    {
        closeBufferedWriter(mGenPosCohortCssWriter);
    }

    private boolean loadSampleCounts()
    {
        if(mSampleDataCache.isSingleSample())
        {
            final String sampleId = mSampleDataCache.SampleIds.get(0);

            // load from VCF, database, sigs file or generic counts files
            if(mConfig.DbAccess != null || !mConfig.PurpleDir.isEmpty())
            {
                List<SomaticVariant> somaticVariants = Lists.newArrayList();

                if(mConfig.DbAccess != null)
                {
                    somaticVariants.addAll(loadSomaticVariants(sampleId, mConfig.DbAccess));
                }
                else
                {
                    String purpleDir = mConfig.getPurpleDataDir(sampleId);
                    final String somaticVcfFile = purpleSomaticVcfFile(purpleDir, sampleId);
                    somaticVariants.addAll(loadSomaticVariants(somaticVcfFile, Lists.newArrayList(SNP)));
                }

                mSampleSnv96Counts = convertSomaticVariantsToSnvCounts(sampleId, somaticVariants, mSampleSnv96CountsIndex);

                AidApobecStatus aidApobecStatus = mIncludeAidApobecGenPos ? AidApobecStatus.ALL : AidApobecStatus.FALSE_ONLY;
                mSampleGenPosCounts = convertSomaticVariantsToPosFrequencies(
                        sampleId, somaticVariants, mSampleGenPosCountsIndex, mPosFrequencies, aidApobecStatus);
            }
            else
            {
                final String snvCountsFile = !mConfig.SampleSnvCountsFile.isEmpty() ?
                        mConfig.SampleSnvCountsFile : mConfig.SampleDataDir + sampleId + ".sig.snv_counts.csv";

                final String snvPosFreqFile = !mConfig.SampleSnvPosFreqFile.isEmpty() ?
                        mConfig.SampleSnvPosFreqFile : mConfig.SampleDataDir + sampleId + ".sig.pos_freq_counts.csv";

                mSampleSnv96Counts = loadSampleCountsFromFile(snvCountsFile, mSampleSnv96CountsIndex);
                mSampleGenPosCounts = loadSampleMatrixData(snvPosFreqFile, mSampleGenPosCountsIndex);
            }

            return mSampleSnv96Counts != null && mSampleGenPosCounts != null;
        }

        if(!mConfig.SampleSnvCountsFile.isEmpty() && !mConfig.SampleSnvPosFreqFile.isEmpty())
        {
            if(mConfig.SampleSnvCountsFile.equals(mConfig.RefSnvCountsFile))
            {
                mSampleSnv96Counts = mRefSampleSnv96Counts;
                mSampleSnv96CountsIndex.putAll(mRefSampleSnv96CountsIndex);
            }
            else
            {
                mSampleSnv96Counts = loadSampleCountsFromFile(mConfig.SampleSnvCountsFile, mSampleSnv96CountsIndex);

                if(mSampleSnv96Counts == null)
                {
                    CUP_LOGGER.error("missing file: {}", mConfig.SampleSnvCountsFile);
                }
            }

            if(mConfig.SampleSnvPosFreqFile.equals(mConfig.RefSnvSamplePosFreqFile))
            {
                mSampleGenPosCounts = mRefSampleGenPosCounts;
                mSampleGenPosCountsIndex.putAll(mRefSampleGenPosCountsIndex);
            }
            else
            {
                mSampleGenPosCounts = loadSampleMatrixData(mConfig.SampleSnvPosFreqFile, mSampleGenPosCountsIndex);

                if(mSampleSnv96Counts == null)
                {
                    CUP_LOGGER.error("missing file: {}", mConfig.SampleSnvPosFreqFile);
                }
            }

            return mSampleSnv96Counts != null && mSampleGenPosCounts != null;
        }

        CUP_LOGGER.error("no sample SNV count source specified");
        return false;
    }

    public void processSample(final SampleData sample, final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        if(!mIsValid || mRefSampleSnv96Counts == null)
            return;

        int snvTotal = mSampleSnvTotals.get(sample.Id);

        addSnv96CssResults(sample, snvTotal, results, similarities);
        addGenomicPositionCssResults(sample, results, similarities);

        mSigContributions.addSigContributionResults(sample, results);

        // add a percentile result
        final Map<String, Double> cancerTypeValues = Maps.newHashMap();

        for(Map.Entry<String, double[]> cancerPercentiles : mRefCancerSnvCountPercentiles.entrySet())
        {
            final String cancerType = cancerPercentiles.getKey();

            if(!isKnownCancerType(cancerType))
                continue;

            double percentile = getPercentile(cancerPercentiles.getValue(), snvTotal, true, UNDEFINED_PERC_MAX_MULTIPLE);
            cancerTypeValues.put(cancerType, percentile);
        }

        SampleResult result = new SampleResult(
                sample.Id, SAMPLE_TRAIT, PERCENTILE, DATA_TYPE_SNV_COUNT, String.valueOf(snvTotal), cancerTypeValues);

        results.add(result);

        int cancerTypeCount = mSampleDataCache.RefCancerSampleData.size();
        int cancerSampleCount = sample.isRefSample() ? mSampleDataCache.getCancerSampleCount(sample.cancerType()) : 0;

        final Map<String,Double> cancerPrevsLow = calcPercentilePrevalence(
                sample, cancerSampleCount, cancerTypeCount, mRefCancerSnvCountPercentiles, snvTotal,  true);

        results.add(new SampleResult(sample.Id, SNV, LIKELIHOOD, DATA_TYPE_SNV_COUNT + "_LOW", String.valueOf(snvTotal), cancerPrevsLow));

        final Map<String,Double> cancerPrevsHigh = calcPercentilePrevalence(
                sample, cancerSampleCount, cancerTypeCount, mRefCancerSnvCountPercentiles, snvTotal, false);

        results.add(new SampleResult(sample.Id, SNV, LIKELIHOOD, DATA_TYPE_SNV_COUNT + "_HIGH", String.valueOf(snvTotal), cancerPrevsHigh));
    }

    private void addSnv96CssResults(
            final SampleData sample, int snvTotal, final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        Integer sampleCountsIndex = mSampleSnv96CountsIndex.get(sample.Id);

        if(sampleCountsIndex == null)
        {
            CUP_LOGGER.info("sample({}) has no SNV data", sample.Id);
            return;
        }

        final double[] sampleCounts = mSampleSnv96Counts.getRow(sampleCountsIndex);

        final List<SampleSimilarity> topMatches = Lists.newArrayList();
        final Map<String,Double> cancerCssTotals = Maps.newHashMap();

        double maxCssScore = 0;

        for(Map.Entry<String,List<SampleData>> refCancerEntry : mSampleDataCache.RefCancerSampleData.entrySet())
        {
            final String refCancerType = refCancerEntry.getKey();

            if(!isKnownCancerType(refCancerType) || !checkIsValidCancerType(sample, refCancerType, cancerCssTotals))
                continue;

            double totalWeightedCss = 0;

            for(SampleData refSample : refCancerEntry.getValue())
            {
                if(refSample.Id.equals(sample.Id))
                    continue;

                Integer refSampleIndex = mRefSampleSnv96CountsIndex.get(refSample.Id);

                if(refSampleIndex == null)
                    continue;

                final double[] otherSampleCounts = mRefSampleSnv96Counts.getRow(refSampleIndex);

                double css = calcCosineSim(sampleCounts, otherSampleCounts);

                if(css < SNV_96_CSS_THRESHOLD)
                    continue;

                if(mWriteSnvSims && mConfig.WriteSimilarities)
                {
                    recordCssSimilarity(
                            topMatches, sample.Id, refSample.Id, css, SNV_96_PAIRWISE.toString(),
                            CSS_SIMILARITY_MAX_MATCHES, CSS_SIMILARITY_CUTOFF);
                }

                maxCssScore = max(css, maxCssScore);

                double cssWeight = pow(mCssExponentSnv, -100 * (1 - css));

                int otherSnvTotal = mSampleSnvTotals.get(refSample.Id);

                double mutLoadWeight = min(otherSnvTotal, snvTotal) / (double)max(otherSnvTotal, snvTotal);

                int cancerTypeCount = mSampleDataCache.getCancerSampleCount(refCancerType);
                double weightedCss = css * cssWeight * mutLoadWeight / sqrt(cancerTypeCount);

                totalWeightedCss += weightedCss;
            }

            cancerCssTotals.put(refCancerType, totalWeightedCss);
        }

        double totalCss = cancerCssTotals.values().stream().mapToDouble(x -> x).sum(); // prior to any conversion

        convertToPercentages(cancerCssTotals);

        if(totalCss > 0)
            applyMaxCssAdjustment(maxCssScore, cancerCssTotals, mMaxCssAdjustFactorSnv);

        results.add(new SampleResult(
                sample.Id, SNV, CLASSIFIER, SNV_96_PAIRWISE.toString(), String.format("%.4g", totalCss), cancerCssTotals));

        // for non-ref cohorts, also report closest matches from amongst these
        if(mWriteSnvSims && mConfig.WriteSimilarities && mSampleDataCache.isMultiSampleNonRef())
        {
            for(Map.Entry<String,Integer> entry : mSampleSnv96CountsIndex.entrySet())
            {
                final String nonRefSampleId = entry.getKey();

                if(nonRefSampleId.equals(sample.Id))
                    continue;

                final double[] otherSampleCounts = mSampleSnv96Counts.getRow(entry.getValue());

                double css = calcCosineSim(sampleCounts, otherSampleCounts);

                if(mConfig.WriteSimilarities)
                {
                    recordCssSimilarity(
                            topMatches, sample.Id, nonRefSampleId, css, SNV_96_PAIRWISE.toString(),
                            CSS_SIMILARITY_MAX_MATCHES, CSS_SIMILARITY_CUTOFF);
                }
            }
        }

        similarities.addAll(topMatches);
    }

    private void addGenomicPositionCssResults(
            final SampleData sample, final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        Integer sampleCountsIndex = mSampleGenPosCountsIndex.get(sample.Id);

        if(sampleCountsIndex == null)
        {
            CUP_LOGGER.debug("sample({}) has no SNV genomic position data", sample.Id);
            return;
        }

        double[] sampleCounts = mSampleGenPosCounts.getRow(sampleCountsIndex);
        double sampleTotal = sumVector(sampleCounts);

        if(mRunPairwiseGenPos)
        {
            addPairwiseGenPosCssResults(sample, sampleCounts, sampleTotal, results, similarities);
        }
        else
        {
            addCohortGenPosCssResults(sample, sampleCounts, sampleTotal, results, similarities);
        }
    }

    private void addCohortGenPosCssResults(
            final SampleData sample, final double[] sampleCounts, double snvTotal,
            final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        // first run CSS against cancer cohorts
        int refCancerCount = mRefCancerGenPosCounts.Rows;
        double maxCssScore = 0;

        final Map<String,Double> cancerCssTotals = Maps.newHashMap();

        for(int i = 0; i < refCancerCount; ++i)
        {
            final String refCancerType = mRefGenPosCancerTypes.get(i);

            if(!isKnownCancerType(refCancerType))
                continue;

            if(!checkIsValidCancerType(sample, refCancerType, cancerCssTotals))
                continue;

            boolean matchesCancerType = sample.cancerType().equals(refCancerType);

            double adjustMultiplier = snvTotal > GEN_POS_MAX_SAMPLE_COUNT ? GEN_POS_MAX_SAMPLE_COUNT / snvTotal : 1;

            final double[] refPosFreqs = sample.isRefSample() && matchesCancerType ?
                    adjustRefCounts(mRefCancerGenPosCounts.getRow(i), sampleCounts, adjustMultiplier) : mRefCancerGenPosCounts.getRow(i);

            double css = calcCosineSim(sampleCounts, refPosFreqs);

            maxCssScore = max(css, maxCssScore);

            writeGenPosCssValues(sample, refCancerType, css);

            cancerCssTotals.put(refCancerType, css);
        }

        double totalWeightedCss = 0;
        for(Map.Entry<String,Double> entry : cancerCssTotals.entrySet())
        {
            String cancerType = entry.getKey();
            double css = entry.getValue();

            double cssWeight;

            if(mGenPosCssExponentCutoff > 0 && mGenPosCssExponentTail > 0)
            {
                cssWeight = pow(mGenPosCssExponent, -100 * min(maxCssScore - css, mGenPosCssExponentCutoff))
                        * pow(mGenPosCssExponentTail, -100 * max(maxCssScore - css, mGenPosCssExponentCutoff));
            }
            else
            {
                cssWeight = pow(mGenPosCssExponent, -100 * (maxCssScore - css));
            }

            double weightedCss = css * cssWeight;

            cancerCssTotals.put(cancerType, weightedCss);
            totalWeightedCss += weightedCss;
        }

        convertToPercentages(cancerCssTotals);

        if(totalWeightedCss > 0)
            applyMaxCssAdjustment(maxCssScore, cancerCssTotals, mMaxCssAdjustFactorGenPos);

        results.add(new SampleResult(
                sample.Id, SNV, CLASSIFIER, GENOMIC_POSITION_COHORT.toString(), String.format("%.4g", totalWeightedCss), cancerCssTotals));
    }


    private void addPairwiseGenPosCssResults(
            final SampleData sample, final double[] sampleCounts, double snvTotal,
            final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        final List<SampleSimilarity> topMatches = Lists.newArrayList();
        final Map<String,Double> cancerCssTotals = Maps.newHashMap();

        double maxCssScore = 0;

        for(Map.Entry<String,List<SampleData>> refCancerEntry : mSampleDataCache.RefCancerSampleData.entrySet())
        {
            final String refCancerType = refCancerEntry.getKey();

            if(!isKnownCancerType(refCancerType) || !checkIsValidCancerType(sample, refCancerType, cancerCssTotals))
                continue;

            double totalWeightedCss = 0;

            for(SampleData refSample : refCancerEntry.getValue())
            {
                if(refSample.Id.equals(sample.Id))
                    continue;

                Integer refSampleIndex = mRefSampleGenPosCountsIndex.get(refSample.Id);

                if(refSampleIndex == null)
                    continue;

                final double[] otherSampleCounts = mRefSampleGenPosCounts.getRow(refSampleIndex);

                double css = calcCosineSim(sampleCounts, otherSampleCounts);

                if(css < GEN_POS_CSS_THRESHOLD)
                    continue;

                if(mWriteGenPosSims && mConfig.WriteSimilarities)
                {
                    recordCssSimilarity(
                            topMatches, sample.Id, refSample.Id, css, GENOMIC_POSITION_COHORT.toString(),
                            20, CSS_SIMILARITY_CUTOFF);

                    /* record scores for all matching cancer-type samples
                    if(sample.cancerType().equals(refCancerType) && topMatches.stream().noneMatch(x -> x.MatchedSampleId.equals(refSample.Id)))
                    {
                        similarities.add(new SampleSimilarity(sample.Id, refSample.Id, GENOMIC_POSITION_SIMILARITY.toString(), css));
                    }
                    */
                }

                if(!isKnownCancerType(refCancerType))
                    continue;

                maxCssScore = max(css, maxCssScore);

                double cssWeight = pow(mGenPosCssExponent, -100 * (1 - css));

                double otherSnvTotal = sumVector(otherSampleCounts);
                double mutLoadWeight = min(otherSnvTotal, snvTotal) / max(otherSnvTotal, snvTotal);

                int cancerTypeCount = mSampleDataCache.getCancerSampleCount(refCancerType);
                double weightedCss = css * cssWeight * mutLoadWeight / sqrt(cancerTypeCount);

                totalWeightedCss += weightedCss;
            }

            cancerCssTotals.put(refCancerType, totalWeightedCss);
        }

        double totalCss = cancerCssTotals.values().stream().mapToDouble(x -> x).sum(); // prior to any conversion

        convertToPercentages(cancerCssTotals);

        if(totalCss > 0)
            applyMaxCssAdjustment(maxCssScore, cancerCssTotals, mMaxCssAdjustFactorSnv);

        results.add(new SampleResult(
                sample.Id, SNV, CLASSIFIER, GENOMIC_POSITION_COHORT.toString(), String.format("%.4g", totalCss), cancerCssTotals));

        similarities.addAll(topMatches);
    }

    private void initialiseOutputFiles()
    {
        try
        {
            final String filename = mConfig.formOutputFilename("GEN_POS_CSS");

            mGenPosCohortCssWriter = createBufferedWriter(filename, false);
            mGenPosCohortCssWriter.write("SampleId,CancerType,RefCancerType,Css");
            mGenPosCohortCssWriter.newLine();

            if(mRefCancerGenPosCounts != null)
            {
                // separately write cancer-type gen pos CSS values
                BufferedWriter writer = createBufferedWriter(mConfig.OutputDir + "CUP.GEN_POS_CANCER_CSS.csv", false);
                writer.write("RefCancerType1,RefCancerType2,Css");
                writer.newLine();

                for(int i = 0; i < mRefGenPosCancerTypes.size() - 1; ++i)
                {
                    final String refCancerType1 = mRefGenPosCancerTypes.get(i);

                    if(!isKnownCancerType(refCancerType1))
                        continue;

                    final double[] refCounts1 = mRefCancerGenPosCounts.getRow(i);

                    for(int j = i + 1; j < mRefGenPosCancerTypes.size(); ++j)
                    {
                        final String refCancerType2 = mRefGenPosCancerTypes.get(j);

                        if(!isKnownCancerType(refCancerType1))
                            continue;

                        final double[] refCounts2 = mRefCancerGenPosCounts.getRow(j);

                        double css = calcCosineSim(refCounts1, refCounts2);

                        writer.write(String.format("%s,%s,%.4f", refCancerType1, refCancerType2, css));
                        writer.newLine();
                    }
                }

                writer.close();
            }
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write gen-pos CSS output: {}", e.toString());
        }
    }

    private synchronized void writeGenPosCssValues(final SampleData sample, final String refCancerType, double css)
    {
        if(mGenPosCohortCssWriter == null)
            return;

        try
        {
            mGenPosCohortCssWriter.write(String.format("%s,%s,%s,%.4f", sample.Id, sample.cancerType(), refCancerType, css));
            mGenPosCohortCssWriter.newLine();
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write sample similarity: {}", e.toString());
        }
    }

    @VisibleForTesting
    public void addRefData(final List<double[]> snvCounts, final List<double[]> posFreqCounts, final Map<String,double[]> cancerPosFreqCounts)
    {
        mRefSampleSnv96Counts = new Matrix(snvCounts.size(), snvCounts.get(0).length);

        for(int i = 0; i < snvCounts.size(); ++i)
        {
            mRefSampleSnv96Counts.setRow(i, snvCounts.get(i));
            mRefSampleSnv96CountsIndex.put(mSampleDataCache.RefSampleDataList.get(i).Id, i);
        }

        mRefCancerGenPosCounts = new Matrix(cancerPosFreqCounts.size(), posFreqCounts.get(0).length);

        int cancerIndex = 0;
        for(Map.Entry<String,double[]> entry : cancerPosFreqCounts.entrySet())
        {
            mRefGenPosCancerTypes.add(entry.getKey());
            mRefCancerGenPosCounts.setRow(cancerIndex, entry.getValue());
            ++cancerIndex;
        }

        mRefSampleGenPosCounts = new Matrix(posFreqCounts.size(), posFreqCounts.get(0).length);

        for(int i = 0; i < posFreqCounts.size(); ++i)
        {
            mRefSampleGenPosCounts.setRow(i, posFreqCounts.get(i));
            mRefSampleGenPosCountsIndex.put(mSampleDataCache.RefSampleDataList.get(i).Id, i);
        }
    }

    public void addSampleData(final List<String> sampleIds, final List<double[]> snvCounts, final List<double[]> posFreqCounts)
    {
        mSampleSnv96Counts = new Matrix(snvCounts.size(), snvCounts.get(0).length);

        for(int i = 0; i < snvCounts.size(); ++i)
        {
            mSampleSnv96Counts.setRow(i, snvCounts.get(i));
            mSampleSnv96CountsIndex.put(sampleIds.get(i), i);
            mSampleSnvTotals.put(sampleIds.get(i), (int)sumVector(snvCounts.get(i)));
        }

        mSampleGenPosCounts = new Matrix(posFreqCounts.size(), posFreqCounts.get(0).length);

        for(int i = 0; i < posFreqCounts.size(); ++i)
        {
            mSampleGenPosCounts.setRow(i, posFreqCounts.get(i));
            mSampleGenPosCountsIndex.put(sampleIds.get(i), i);
        }
    }

}
