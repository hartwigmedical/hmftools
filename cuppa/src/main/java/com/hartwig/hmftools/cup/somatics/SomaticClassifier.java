package com.hartwig.hmftools.cup.somatics;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;
import static java.lang.Math.sqrt;

import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;
import static com.hartwig.hmftools.common.stats.CosineSimilarity.calcCosineSim;
import static com.hartwig.hmftools.common.stats.Percentiles.getPercentile;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.VariantType.SNP;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CategoryType.CLASSIFIER;
import static com.hartwig.hmftools.cup.common.CategoryType.SAMPLE_TRAIT;
import static com.hartwig.hmftools.cup.common.CategoryType.SNV;
import static com.hartwig.hmftools.cup.common.ClassifierType.GENOMIC_POSITION_SIMILARITY;
import static com.hartwig.hmftools.cup.common.ClassifierType.SNV_96_PAIRWISE_SIMILARITY;
import static com.hartwig.hmftools.cup.common.CupCalcs.adjustRefCounts;
import static com.hartwig.hmftools.cup.common.CupCalcs.calcPercentilePrevalence;
import static com.hartwig.hmftools.cup.common.CupCalcs.convertToPercentages;
import static com.hartwig.hmftools.cup.common.CupConstants.AID_APOBEC_TRINUCLEOTIDE_CONTEXTS;
import static com.hartwig.hmftools.cup.common.CupConstants.CSS_SIMILARITY_CUTOFF;
import static com.hartwig.hmftools.cup.common.CupConstants.CSS_SIMILARITY_MAX_MATCHES;
import static com.hartwig.hmftools.cup.common.CupConstants.POS_FREQ_BUCKET_SIZE;
import static com.hartwig.hmftools.cup.common.CupConstants.POS_FREQ_MAX_SAMPLE_COUNT;
import static com.hartwig.hmftools.cup.common.CupConstants.SNV_CSS_DIFF_EXPONENT;
import static com.hartwig.hmftools.cup.common.CupConstants.SNV_CSS_THRESHOLD;
import static com.hartwig.hmftools.cup.common.CupConstants.SNV_POS_FREQ_CSS_THRESHOLD;
import static com.hartwig.hmftools.cup.common.CupConstants.SNV_POS_FREQ_DIFF_EXPONENT;
import static com.hartwig.hmftools.cup.common.CupConstants.UNDEFINED_PERC_MAX_MULTIPLE;
import static com.hartwig.hmftools.cup.common.ResultType.LIKELIHOOD;
import static com.hartwig.hmftools.cup.common.ResultType.PERCENTILE;
import static com.hartwig.hmftools.cup.common.SampleData.isKnownCancerType;
import static com.hartwig.hmftools.cup.common.SampleResult.checkIsValidCancerType;
import static com.hartwig.hmftools.cup.common.SampleSimilarity.recordCssSimilarity;
import static com.hartwig.hmftools.cup.somatics.CopyNumberProfile.normaliseGenPosCountsByCopyNumber;
import static com.hartwig.hmftools.cup.somatics.GenomicPositions.convertSomaticVariantsToPosFrequencies;
import static com.hartwig.hmftools.cup.somatics.SomaticDataLoader.loadRefSampleCounts;
import static com.hartwig.hmftools.cup.somatics.SomaticDataLoader.loadRefSignaturePercentileData;
import static com.hartwig.hmftools.cup.somatics.SomaticDataLoader.loadSampleCountsFromFile;
import static com.hartwig.hmftools.cup.somatics.SomaticDataLoader.loadSampleMatrixData;
import static com.hartwig.hmftools.cup.somatics.SomaticDataLoader.loadSomaticVariants;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.EXCLUDE_SNV_96_AID_APOBEC;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.EXCLUDE_SNV_96_AID_APOBEC_DESC;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.INCLUDE_AID_APOBEC_SIG;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.INCLUDE_AID_APOBEC_SIG_DESC;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.NORMALISE_COPY_NUMBER;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.NORMALISE_COPY_NUMBER_DESC;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.SPLIT_AID_APOBEC;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.SPLIT_AID_APOBEC_DESC;
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
import com.hartwig.hmftools.common.sigs.SnvSigUtils;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.cup.CuppaConfig;
import com.hartwig.hmftools.cup.common.CategoryType;
import com.hartwig.hmftools.cup.common.CuppaClassifier;
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

    private Matrix mRefSampleCounts;
    private final List<String> mRefSampleNames;
    private final Map<String,double[]> mRefCancerSnvCountPercentiles;

    private Matrix mRefCancerSnvPosFrequencies;
    private final List<String> mRefSnvPosFreqCancerTypes;

    private Matrix mRefSamplePosFrequencies;
    private final Map<String,Integer> mRefSamplePosFreqIndex;

    private Matrix mSampleSnvCounts;
    private final Map<String,Integer> mSampleSnvCountsIndex;

    private Matrix mSamplePosFrequencies;
    private final Map<String,Integer> mSamplePosFreqIndex;

    private final SigContributions mSigContributions;
    private final PositionFrequencies mPosFrequencies;

    private boolean mIsValid;
    private BufferedWriter mCssWriter;

    private final boolean mSplitAidApobecGenPos;
    private final boolean mExcludeAidApobecSnv96;

    private final boolean mApplyCopyNumber; // to genomic positions
    private Matrix mRefSampleCopyNumberProfiles;
    private Matrix mCopyNumberProfile; // same dimensions as genomic position

    private final double mMaxCssAdjustFactorSnv;
    private final double mMaxCssAdjustFactorGenPos;
    private final double mCssExponentSnv;
    private final double mCssExponentGenPos;
    private final boolean mWriteGenPosSims;
    private final boolean mWriteSnvSims;

    private final List<Integer> mAidApobecSnv96Buckets;

    public static final String MAX_CSS_ADJUST_FACTOR_SNV = "css_max_factor_snv";
    public static final String MAX_CSS_ADJUST_FACTOR_GEN_POS = "css_max_factor_gen_pos";

    public static final String CSS_EXPONENT_SNV = "css_exponent_snv";
    public static final String CSS_EXPONENT_GEN_POS = "css_exponent_gen_pos";

    public static final String SNV_POS_FREQ_POS_SIZE = "pos_freq_bucket_size";

    public static final String WRITE_GEN_POS_CSS = "write_gen_pos_css";
    public static final String WRITE_GEN_POS_SIMILARITIES = "write_gen_pos_sims";
    public static final String WRITE_SNV_SIMILARITIES = "write_snv_sims";

    private static final int GEN_POS_CSS_SIMILARITY_MAX_MATCHES = 100;

    public SomaticClassifier(final CuppaConfig config, final SampleDataCache sampleDataCache, final CommandLine cmd)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mSampleSnvCounts = null;
        mRefCancerSnvPosFrequencies = null;
        mRefSamplePosFrequencies = null;
        mRefSampleCopyNumberProfiles = null;
        mSampleSnvCountsIndex = Maps.newHashMap();
        mSamplePosFreqIndex = Maps.newHashMap();

        mRefSampleCounts = null;
        mRefSampleNames = Lists.newArrayList();
        mRefCancerSnvCountPercentiles = Maps.newHashMap();
        mRefSnvPosFreqCancerTypes = Lists.newArrayList();
        mRefSamplePosFreqIndex = Maps.newHashMap();

        mSplitAidApobecGenPos = cmd != null ? cmd.hasOption(SPLIT_AID_APOBEC) : false;
        mExcludeAidApobecSnv96 = cmd != null ? cmd.hasOption(EXCLUDE_SNV_96_AID_APOBEC) : false;
        mApplyCopyNumber = cmd != null ? cmd.hasOption(NORMALISE_COPY_NUMBER) : false;
        mAidApobecSnv96Buckets = Lists.newArrayList();

        mCssExponentSnv = cmd != null ? Double.parseDouble(cmd.getOptionValue(CSS_EXPONENT_SNV, "8")) : SNV_CSS_DIFF_EXPONENT;
        mCssExponentGenPos = cmd != null ? Double.parseDouble(cmd.getOptionValue(CSS_EXPONENT_GEN_POS, "10")) : SNV_POS_FREQ_DIFF_EXPONENT;
        mMaxCssAdjustFactorSnv = cmd != null ? Double.parseDouble(cmd.getOptionValue(MAX_CSS_ADJUST_FACTOR_SNV, "0")) : 0;
        mMaxCssAdjustFactorGenPos = cmd != null ? Double.parseDouble(cmd.getOptionValue(MAX_CSS_ADJUST_FACTOR_GEN_POS, "0")) : 0;
        mWriteSnvSims = cmd != null ? cmd.hasOption(WRITE_SNV_SIMILARITIES) : false;
        mWriteGenPosSims = cmd != null ? cmd.hasOption(WRITE_GEN_POS_SIMILARITIES) : false;

        mIsValid = true;
        mCssWriter = null;

        mSigContributions = new SigContributions(mConfig, mSampleDataCache);

        int posFreqBucketSize = cmd != null && cmd.hasOption(SNV_POS_FREQ_POS_SIZE) ?
                Integer.parseInt(cmd.getOptionValue(SNV_POS_FREQ_POS_SIZE)) : POS_FREQ_BUCKET_SIZE;

        mPosFrequencies = new PositionFrequencies(posFreqBucketSize, POS_FREQ_MAX_SAMPLE_COUNT);

        if(mConfig.RefSnvCountsFile.isEmpty() && mConfig.RefSigContributionFile.isEmpty() && mConfig.RefSnvCancerPosFreqFile.isEmpty())
            return;

        loadRefSignaturePercentileData(
                mConfig.RefSigContributionFile, mSigContributions.getRefCancerSigContribPercentiles(), mRefCancerSnvCountPercentiles);

        mRefSampleCounts = loadRefSampleCounts(mConfig.RefSnvCountsFile, mRefSampleNames, Lists.newArrayList("BucketName"));

        mRefCancerSnvPosFrequencies = loadRefSampleCounts(mConfig.RefSnvCancerPosFreqFile, mRefSnvPosFreqCancerTypes, Lists.newArrayList());
        mRefSamplePosFrequencies = loadSampleMatrixData(mConfig.RefSnvSamplePosFreqFile, mRefSamplePosFreqIndex);

        if(mApplyCopyNumber)
        {
            // uses gen-pos sample map for CN profile
            mRefSampleCopyNumberProfiles = loadSampleMatrixData(mConfig.RefCopyNumberProfileFile, Maps.newHashMap());

            if(!mConfig.SampleDataFile.equals(mConfig.RefSampleDataFile))
            {
                CUP_LOGGER.error("only ref-sample analysis support for CN normalisation");
                mIsValid = false;
            }
        }
        else
        {
            mRefSampleCopyNumberProfiles = null;
        }

        if(mRefSampleCounts == null || mRefCancerSnvPosFrequencies == null || (mApplyCopyNumber && mRefSampleCopyNumberProfiles == null))
        {
            CUP_LOGGER.error("invalid somatic matrix data: SNV-96() GenPos({}) CopyNumber({})",
                    mRefSampleCounts != null, mRefCancerSnvPosFrequencies != null,
                    (!mApplyCopyNumber || mRefSampleCopyNumberProfiles != null));
            mIsValid = false;
        }

        mIsValid &= loadSampleCounts();
        mIsValid &= mSigContributions.loadSigContributions(mSampleSnvCounts);

        if(mApplyCopyNumber && !mConfig.runClassifier(SAMPLE_TRAIT))
        {
            CUP_LOGGER.error("gennomic-position copy number normalisation requires sample traits for purity");
            mIsValid = false;
        }

        excludeAidApobecBuckets();

        if(cmd.hasOption(WRITE_GEN_POS_CSS))
        {
            initialiseOutputFiles();
        }
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(SPLIT_AID_APOBEC, false, SPLIT_AID_APOBEC_DESC);
        options.addOption(NORMALISE_COPY_NUMBER, false, NORMALISE_COPY_NUMBER_DESC);
        options.addOption(EXCLUDE_SNV_96_AID_APOBEC, false, EXCLUDE_SNV_96_AID_APOBEC_DESC);
        options.addOption(INCLUDE_AID_APOBEC_SIG, false, INCLUDE_AID_APOBEC_SIG_DESC);
        options.addOption(MAX_CSS_ADJUST_FACTOR_SNV, true, "Max CSS adustment factor for SNV 96");
        options.addOption(MAX_CSS_ADJUST_FACTOR_GEN_POS, true, "Max CSS adustment factor for genomic pos frequency");
        options.addOption(CSS_EXPONENT_SNV, true, "Max CSS adustment factor for SNV 96");
        options.addOption(CSS_EXPONENT_GEN_POS, true, "Max CSS adustment factor for SNV 96");
        options.addOption(SNV_POS_FREQ_POS_SIZE, true, "Genomic position bucket size (default: 20000)");
        options.addOption(WRITE_GEN_POS_CSS, false, "Write gen-pos CSS to file");
    }

    public CategoryType categoryType() { return SNV; }
    public boolean isValid() { return mIsValid; }

    public void close()
    {
        closeBufferedWriter(mCssWriter);
    }

    private boolean loadSampleCounts()
    {
        if(mSampleDataCache.isSingleSample())
        {
            final String sampleId = mSampleDataCache.SampleIds.get(0);

            // load from VCF, database, sigs file or generic counts files
            if(mConfig.DbAccess != null || !mConfig.SampleSomaticVcf.isEmpty())
            {
                List<SomaticVariant> somaticVariants = Lists.newArrayList();

                if(mConfig.DbAccess != null)
                {
                    somaticVariants.addAll(loadSomaticVariants(sampleId, mConfig.DbAccess));
                }
                else
                {
                    somaticVariants.addAll(loadSomaticVariants(mConfig.SampleSomaticVcf, Lists.newArrayList(SNP)));
                }

                mSampleSnvCounts = convertSomaticVariantsToSnvCounts(sampleId, somaticVariants, mSampleSnvCountsIndex);

                AidApobecStatus aidApobecStatus = mSplitAidApobecGenPos ? AidApobecStatus.FALSE_ONLY : AidApobecStatus.ALL;
                mSamplePosFrequencies = convertSomaticVariantsToPosFrequencies(
                        sampleId, somaticVariants, mSamplePosFreqIndex, mPosFrequencies, aidApobecStatus);
            }
            else
            {
                final String snvCountsFile = !mConfig.SampleSnvCountsFile.isEmpty() ?
                        mConfig.SampleSnvCountsFile : mConfig.SampleDataDir + sampleId + ".sig.snv_counts.csv";

                final String snvPosFreqFile = !mConfig.SampleSnvPosFreqFile.isEmpty() ?
                        mConfig.SampleSnvPosFreqFile : mConfig.SampleDataDir + sampleId + ".sig.pos_freq_counts.csv";

                mSampleSnvCounts = loadSampleCountsFromFile(snvCountsFile, mSampleSnvCountsIndex);
                mSamplePosFrequencies = loadSampleMatrixData(snvPosFreqFile, mSamplePosFreqIndex);
            }

            return mSampleSnvCounts != null && mSamplePosFrequencies != null;
        }

        if(!mConfig.SampleSnvCountsFile.isEmpty() && !mConfig.SampleSnvPosFreqFile.isEmpty())
        {
            if(mConfig.SampleSnvCountsFile.equals(mConfig.RefSnvCountsFile))
            {
                mSampleSnvCounts = mRefSampleCounts;

                for(int i = 0; i < mRefSampleNames.size(); ++i)
                {
                    mSampleSnvCountsIndex.put(mRefSampleNames.get(i), i);
                }
            }
            else
            {
                mSampleSnvCounts = loadSampleCountsFromFile(mConfig.SampleSnvCountsFile, mSampleSnvCountsIndex);
            }

            if(mConfig.SampleSnvPosFreqFile.equals(mConfig.RefSnvSamplePosFreqFile))
            {
                mSamplePosFrequencies = mRefSamplePosFrequencies;
                mSamplePosFreqIndex.putAll(mRefSamplePosFreqIndex);
            }
            else
            {
                mSamplePosFrequencies = loadSampleMatrixData(mConfig.SampleSnvPosFreqFile, mSamplePosFreqIndex);
            }

            return mSampleSnvCounts != null && mSamplePosFrequencies != null;
        }

        CUP_LOGGER.error("no sample SNV count source specified");
        return false;
    }

    public void processSample(final SampleData sample, final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        if(!mIsValid || mRefSampleCounts == null)
            return;

        Integer sampleCountsIndex = mSampleSnvCountsIndex.get(sample.Id);

        if(sampleCountsIndex == null)
        {
            CUP_LOGGER.info("sample({}) has no SNV data", sample.Id);
            return;
        }

        final double[] sampleCounts = mSampleSnvCounts.getCol(sampleCountsIndex);
        int snvTotal = (int)sumVector(sampleCounts);

        addSnv96CssResults(sample, sampleCounts, snvTotal, results, similarities);
        addPosFreqCssResults(sample, results, similarities);

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
                sample.Id, SAMPLE_TRAIT, PERCENTILE, "SNV_COUNT", String.valueOf(snvTotal), cancerTypeValues);

        results.add(result);

        int cancerTypeCount = mSampleDataCache.RefCancerSampleData.size();
        int cancerSampleCount = sample.isRefSample() ? mSampleDataCache.getCancerSampleCount(sample.cancerType()) : 0;

        final Map<String,Double> cancerPrevsLow = calcPercentilePrevalence(
                sample, cancerSampleCount, cancerTypeCount, mRefCancerSnvCountPercentiles, snvTotal,  true);

        results.add(new SampleResult(sample.Id, SNV, LIKELIHOOD, "SNV_COUNT_LOW", String.valueOf(snvTotal), cancerPrevsLow));

        final Map<String,Double> cancerPrevsHigh = calcPercentilePrevalence(
                sample, cancerSampleCount, cancerTypeCount, mRefCancerSnvCountPercentiles, snvTotal, false);

        results.add(new SampleResult(sample.Id, SNV, LIKELIHOOD, "SNV_COUNT_HIGH", String.valueOf(snvTotal), cancerPrevsHigh));
    }

    private void addSnv96CssResults(
            final SampleData sample, final double[] sampleCounts, int snvTotal,
            final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        int refSampleCount = mRefSampleCounts.Cols;

        final List<SampleSimilarity> topMatches = Lists.newArrayList();
        final Map<String,Double> cancerCssTotals = Maps.newHashMap();

        if(mExcludeAidApobecSnv96)
        {
            for(Integer bucketIndex : mAidApobecSnv96Buckets)
            {
                sampleCounts[bucketIndex] = 0;
            }
        }

        double maxCssScore = 0;

        for(int s = 0; s < refSampleCount; ++s)
        {
            final String refSampleId = mRefSampleNames.get(s);

            if(!mSampleDataCache.hasRefSample(refSampleId))
                continue;

            if(refSampleId.equals(sample.Id))
                continue;

            final String refCancerType = mSampleDataCache.RefSampleCancerTypeMap.get(refSampleId);

            if(!checkIsValidCancerType(sample, refCancerType, cancerCssTotals))
                continue;

            final double[] otherSampleCounts = mRefSampleCounts.getCol(s);

            double css = calcCosineSim(sampleCounts, otherSampleCounts);

            if(css < SNV_CSS_THRESHOLD)
                continue;

            if(mWriteSnvSims && mConfig.WriteSimilarities)
            {
                recordCssSimilarity(
                        topMatches, sample.Id, refSampleId, css, SNV_96_PAIRWISE_SIMILARITY.toString(),
                        CSS_SIMILARITY_MAX_MATCHES, CSS_SIMILARITY_CUTOFF);
            }

            if(!isKnownCancerType(refCancerType))
                continue;

            maxCssScore = max(css, maxCssScore);

            double cssWeight = pow(mCssExponentSnv, -100 * (1 - css));

            double otherSnvTotal = sumVector(otherSampleCounts);
            double mutLoadWeight = min(otherSnvTotal, snvTotal) / max(otherSnvTotal, snvTotal);

            int cancerTypeCount = mSampleDataCache.getCancerSampleCount(refCancerType);
            double weightedCss = css * cssWeight * mutLoadWeight / sqrt(cancerTypeCount);

            Double total = cancerCssTotals.get(refCancerType);

            if(total == null)
                cancerCssTotals.put(refCancerType, weightedCss);
            else
                cancerCssTotals.put(refCancerType, total + weightedCss);
        }

        double totalCss = cancerCssTotals.values().stream().mapToDouble(x -> x).sum(); // prior to any conversion

        convertToPercentages(cancerCssTotals);

        if(totalCss > 0)
            applyMaxCssAdjustment(maxCssScore, cancerCssTotals, mMaxCssAdjustFactorSnv);

        results.add(new SampleResult(
                sample.Id, CLASSIFIER, LIKELIHOOD, SNV_96_PAIRWISE_SIMILARITY.toString(), String.format("%.4g", totalCss), cancerCssTotals));

        // for non-ref cohorts, also report closest matches from amongst these
        if(mWriteSnvSims && mConfig.WriteSimilarities && mSampleDataCache.isMultiSampleNonRef())
        {
            for(Map.Entry<String,Integer> entry : mSampleSnvCountsIndex.entrySet())
            {
                final String nonRefSampleId = entry.getKey();

                if(nonRefSampleId.equals(sample.Id))
                    continue;

                final double[] otherSampleCounts = mSampleSnvCounts.getCol(entry.getValue());

                double css = calcCosineSim(sampleCounts, otherSampleCounts);

                if(mConfig.WriteSimilarities)
                {
                    recordCssSimilarity(
                            topMatches, sample.Id, nonRefSampleId, css, SNV_96_PAIRWISE_SIMILARITY.toString(),
                            CSS_SIMILARITY_MAX_MATCHES, CSS_SIMILARITY_CUTOFF);
                }
            }
        }

        similarities.addAll(topMatches);
    }

    private void addPosFreqCssResults(
            final SampleData sample, final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        Integer sampleCountsIndex = mSamplePosFreqIndex.get(sample.Id);

        if(sampleCountsIndex == null)
        {
            CUP_LOGGER.debug("sample({}) has no SNV pos-freq data", sample.Id);
            return;
        }

        double[] sampleCounts = mSamplePosFrequencies.getCol(sampleCountsIndex);

        if(mApplyCopyNumber)
        {
            double samplePloidy = mSampleDataCache.RefSampleTraitsData.get(sample.Id).Ploidy;
            final double[] sampleCnProfile = mRefSampleCopyNumberProfiles.getCol(sampleCountsIndex);
            sampleCounts = normaliseGenPosCountsByCopyNumber(samplePloidy, sampleCounts, sampleCnProfile);
        }

        double sampleTotal = sumVector(sampleCounts);

        // first run CSS against cancer cohorts
        int refCancerCount = mRefCancerSnvPosFrequencies.Cols;
        double maxCssScore = 0;

        final Map<String,Double> cancerCssTotals = Maps.newHashMap();

        for(int i = 0; i < refCancerCount; ++i)
        {
            final String refCancerType = mRefSnvPosFreqCancerTypes.get(i);

            if(!isKnownCancerType(refCancerType))
                continue;

            if(!checkIsValidCancerType(sample, refCancerType, cancerCssTotals))
                continue;

            boolean matchesCancerType = sample.cancerType().equals(refCancerType);

            double adjustMultiplier = sampleTotal > POS_FREQ_MAX_SAMPLE_COUNT ? POS_FREQ_MAX_SAMPLE_COUNT / sampleTotal : 1;

            final double[] refPosFreqs = sample.isRefSample() && matchesCancerType ?
                    adjustRefCounts(mRefCancerSnvPosFrequencies.getCol(i), sampleCounts, adjustMultiplier) : mRefCancerSnvPosFrequencies.getCol(i);

            double css = calcCosineSim(sampleCounts, refPosFreqs);

            maxCssScore = max(css, maxCssScore);

            writeGenPosCssValues(sample, refCancerType, css);

            if(css < SNV_POS_FREQ_CSS_THRESHOLD)
                continue;

            double cssWeight = pow(mCssExponentGenPos, -100 * (1 - css));

            double weightedCss = css * cssWeight;

            Double total = cancerCssTotals.get(refCancerType);

            if(total == null)
                cancerCssTotals.put(refCancerType, weightedCss);
            else
                cancerCssTotals.put(refCancerType, total + weightedCss);
        }

        double totalCss = cancerCssTotals.values().stream().mapToDouble(x -> x).sum();

        convertToPercentages(cancerCssTotals);

        if(totalCss > 0)
            applyMaxCssAdjustment(maxCssScore, cancerCssTotals, mMaxCssAdjustFactorGenPos);

        results.add(new SampleResult(
                sample.Id, CLASSIFIER, LIKELIHOOD, GENOMIC_POSITION_SIMILARITY.toString(), String.format("%.4g", totalCss), cancerCssTotals));

        // then run pairwise analysis if similarities are being analysed
        if(mWriteGenPosSims)
            addSnvPosSimilarities(sample, sampleCounts, similarities);
    }

    private void initialiseOutputFiles()
    {
        try
        {
            final String filename = mConfig.OutputDir + "CUP.GEN_POS_CSS.csv";

            mCssWriter = createBufferedWriter(filename, false);
            mCssWriter.write("SampleId,CancerType,RefCancerType,Css");
            mCssWriter.newLine();

            if(mRefCancerSnvPosFrequencies != null)
            {
                // separately write cancer-type gen pos CSS values
                BufferedWriter writer = createBufferedWriter(mConfig.OutputDir + "CUP.GEN_POS_CANCER_CSS.csv", false);
                writer.write("RefCancerType1,RefCancerType2,Css");
                writer.newLine();

                for(int i = 0; i < mRefSnvPosFreqCancerTypes.size() - 1; ++i)
                {
                    final String refCancerType1 = mRefSnvPosFreqCancerTypes.get(i);

                    if(!isKnownCancerType(refCancerType1))
                        continue;

                    final double[] refCounts1 = mRefCancerSnvPosFrequencies.getCol(i);

                    for(int j = i + 1; j < mRefSnvPosFreqCancerTypes.size(); ++j)
                    {
                        final String refCancerType2 = mRefSnvPosFreqCancerTypes.get(j);

                        if(!isKnownCancerType(refCancerType1))
                            continue;

                        final double[] refCounts2 = mRefCancerSnvPosFrequencies.getCol(j);

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

    private void writeGenPosCssValues(final SampleData sample, final String refCancerType, double css)
    {
        if(mCssWriter == null)
            return;

        try
        {
            mCssWriter.write(String.format("%s,%s,%s,%.4f", sample.Id, sample.cancerType(), refCancerType, css));
            mCssWriter.newLine();
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write sample similarity: {}", e.toString());
        }
    }

    private void addSnvPosSimilarities(final SampleData sample, final double[] sampleCounts, final List<SampleSimilarity> similarities)
    {
        // not currently used
        if(mRefSamplePosFrequencies == null || !mConfig.WriteSimilarities)
            return;

        final List<SampleSimilarity> topMatches = Lists.newArrayList();

        for(Map.Entry<String,Integer> entry : mRefSamplePosFreqIndex.entrySet())
        {
            final String refSampleId = entry.getKey();

            if(refSampleId.equals(sample.Id))
                continue;

            final String refCancerType = mSampleDataCache.RefSampleCancerTypeMap.get(refSampleId);

            if(refCancerType == null || !sample.isCandidateCancerType(refCancerType))
                continue;

            final double[] otherSampleCounts = mRefSamplePosFrequencies.getCol(entry.getValue());

            double css = calcCosineSim(sampleCounts, otherSampleCounts);

            if(css < SNV_CSS_THRESHOLD)
                continue;

            recordCssSimilarity(
                    topMatches, sample.Id, refSampleId, css, GENOMIC_POSITION_SIMILARITY.toString(),
                    GEN_POS_CSS_SIMILARITY_MAX_MATCHES, CSS_SIMILARITY_CUTOFF);
        }

        similarities.addAll(topMatches);
    }

    private void excludeAidApobecBuckets()
    {
        if(!mExcludeAidApobecSnv96)
            return;

        Map<String,Integer> bucketNameIndexMap = Maps.newHashMap();
        SnvSigUtils.populateBucketMap(bucketNameIndexMap);
        for(String bucketName : AID_APOBEC_TRINUCLEOTIDE_CONTEXTS)
        {
            int bucketIndex = bucketNameIndexMap.get(bucketName);
            mAidApobecSnv96Buckets.add(bucketIndex);

            for(int i = 0; i < mRefSampleCounts.Cols; ++i)
            {
                mRefSampleCounts.set(bucketIndex, i, 0);
            }
        }
    }

    @VisibleForTesting
    public void addRefData(final List<double[]> snvCounts, final List<double[]> posFreqCounts, final Map<String,double[]> cancerPosFreqCounts)
    {
        mRefSampleCounts = new Matrix(snvCounts.get(0).length, snvCounts.size());

        for(int i = 0; i < snvCounts.size(); ++i)
        {
            mRefSampleCounts.setCol(i, snvCounts.get(i));
            mRefSampleNames.add(mSampleDataCache.RefSampleDataList.get(i).Id);
        }

        mRefCancerSnvPosFrequencies = new Matrix(posFreqCounts.get(0).length, cancerPosFreqCounts.size());

        int cancerIndex = 0;
        for(Map.Entry<String,double[]> entry : cancerPosFreqCounts.entrySet())
        {
            mRefSnvPosFreqCancerTypes.add(entry.getKey());
            mRefCancerSnvPosFrequencies.setCol(cancerIndex, entry.getValue());
            ++cancerIndex;
        }

        mRefSamplePosFrequencies = new Matrix(posFreqCounts.get(0).length, posFreqCounts.size());

        for(int i = 0; i < posFreqCounts.size(); ++i)
        {
            mRefSamplePosFrequencies.setCol(i, posFreqCounts.get(i));
            mRefSamplePosFreqIndex.put(mSampleDataCache.RefSampleDataList.get(i).Id, i);
        }
    }

    public void addSampleData(final List<String> sampleIds, final List<double[]> snvCounts, final List<double[]> posFreqCounts)
    {
        mSampleSnvCounts = new Matrix(snvCounts.get(0).length, snvCounts.size());

        for(int i = 0; i < snvCounts.size(); ++i)
        {
            mSampleSnvCounts.setCol(i, snvCounts.get(i));
            mSampleSnvCountsIndex.put(sampleIds.get(i), i);
        }

        mSamplePosFrequencies = new Matrix(posFreqCounts.get(0).length, posFreqCounts.size());

        for(int i = 0; i < posFreqCounts.size(); ++i)
        {
            mSamplePosFrequencies.setCol(i, posFreqCounts.get(i));
            mSamplePosFreqIndex.put(sampleIds.get(i), i);
        }
    }

}
