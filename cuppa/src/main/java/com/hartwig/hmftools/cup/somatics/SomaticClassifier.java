package com.hartwig.hmftools.cup.somatics;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

import static com.hartwig.hmftools.common.sigs.SnvSigUtils.populateBucketMap;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;
import static com.hartwig.hmftools.common.stats.CosineSimilarity.calcCosineSim;
import static com.hartwig.hmftools.common.stats.Percentiles.getPercentile;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.VariantType.SNP;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.common.cuppa.CategoryType.SAMPLE_TRAIT;
import static com.hartwig.hmftools.common.cuppa.CategoryType.SNV;
import static com.hartwig.hmftools.common.cuppa.ClassifierType.GENOMIC_POSITION_COHORT;
import static com.hartwig.hmftools.common.cuppa.ClassifierType.SNV_96_PAIRWISE;
import static com.hartwig.hmftools.cup.CuppaRefFiles.purpleSomaticVcfFile;
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
import static com.hartwig.hmftools.common.cuppa.ResultType.CLASSIFIER;
import static com.hartwig.hmftools.common.cuppa.ResultType.LIKELIHOOD;
import static com.hartwig.hmftools.common.cuppa.ResultType.PERCENTILE;
import static com.hartwig.hmftools.cup.common.SampleData.isKnownCancerType;
import static com.hartwig.hmftools.cup.common.SampleResult.checkIsValidCancerType;
import static com.hartwig.hmftools.cup.common.SampleSimilarity.recordCssSimilarity;
import static com.hartwig.hmftools.cup.somatics.GenomicPositions.extractPositionFrequencyCounts;
import static com.hartwig.hmftools.cup.somatics.SomaticDataLoader.loadRefSampleCounts;
import static com.hartwig.hmftools.cup.somatics.SomaticDataLoader.loadRefSignaturePercentileData;
import static com.hartwig.hmftools.cup.somatics.SomaticDataLoader.loadSampleCountsFromFile;
import static com.hartwig.hmftools.cup.somatics.SomaticDataLoader.loadSampleMatrixData;
import static com.hartwig.hmftools.cup.somatics.SomaticDataLoader.loadSomaticVariants;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.GEN_POS_BUCKET_SIZE_CFG;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.GEN_POS_BUCKET_SIZE_DESC;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.GEN_POS_MAX_SAMPLE_COUNT_CFG;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.GEN_POS_MAX_SAMPLE_COUNT_DESC;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.INCLUDE_AID_APOBEC_SIG;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.INCLUDE_AID_APOBEC_SIG_DESC;
import static com.hartwig.hmftools.cup.somatics.TrinucleotideCounts.extractTrinucleotideCounts;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sigs.PositionFrequencies;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.cup.CuppaConfig;
import com.hartwig.hmftools.common.cuppa.CategoryType;
import com.hartwig.hmftools.cup.common.CuppaClassifier;
import com.hartwig.hmftools.cup.common.NoiseRefCache;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;
import com.hartwig.hmftools.cup.common.SampleSimilarity;

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

    private BufferedWriter mGenPosCohortCssWriter;

    private final double mCssExponentSnv;
    private final double mGenPosCssExponent;
    private final double mGenPosCssExponentTail;
    private final double mGenPosCssExponentCutoff;
    private final boolean mWriteGenPosSims;
    private final boolean mWriteSnvSims;

    // config
    public static final String CSS_EXPONENT_SNV = "snv_css_exponent";
    public static final String CSS_EXPONENT_GEN_POS = "gen_pos_css_exponent";
    public static final String CSS_EXPONENT_GEN_POS_TAIL = "gen_pos_css_exponent_tail";
    public static final String CSS_EXPONENT_GEN_POS_CUTOFF = "gen_pos_css_exponent_cutoff";

    public static final String WRITE_GEN_POS_CSS = "write_gen_pos_css";
    public static final String WRITE_GEN_POS_SIMILARITIES = "write_gen_pos_sims";
    public static final String WRITE_SNV_SIMILARITIES = "write_snv_sims";

    public SomaticClassifier(final CuppaConfig config, final SampleDataCache sampleDataCache, final ConfigBuilder configBuilder)
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

        mCssExponentSnv = configBuilder.getDecimal(CSS_EXPONENT_SNV);
        mGenPosCssExponent = configBuilder.getDecimal(CSS_EXPONENT_GEN_POS);
        mGenPosCssExponentTail = configBuilder.getDecimal(CSS_EXPONENT_GEN_POS_TAIL);
        mGenPosCssExponentCutoff = configBuilder.getDecimal(CSS_EXPONENT_GEN_POS_CUTOFF);

        mWriteSnvSims = configBuilder.hasFlag(WRITE_SNV_SIMILARITIES);
        mWriteGenPosSims = configBuilder.hasFlag(WRITE_GEN_POS_SIMILARITIES);

        mGenPosCohortCssWriter = null;

        mSigContributions = new SigContributions(mConfig, mSampleDataCache);

        int genPosBucketSize = configBuilder.getInteger(GEN_POS_BUCKET_SIZE_CFG);
        int genPosMaxSampleCount = configBuilder.getInteger(GEN_POS_MAX_SAMPLE_COUNT_CFG);

        mPosFrequencies = new PositionFrequencies(mConfig.RefGenVersion, genPosBucketSize, genPosMaxSampleCount);

        if(configBuilder.hasFlag(WRITE_GEN_POS_CSS))
        {
            initialiseOutputFiles();
        }
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addFlag(INCLUDE_AID_APOBEC_SIG, INCLUDE_AID_APOBEC_SIG_DESC);
        configBuilder.addDecimal(CSS_EXPONENT_SNV, "SNV 96 CSS exponent for tail", SNV_96_CSS_DIFF_EXPONENT);
        configBuilder.addDecimal(CSS_EXPONENT_GEN_POS, "Genomic position CSS exponent", GEN_POS_CSS_EXPONENT);
        configBuilder.addDecimal(CSS_EXPONENT_GEN_POS_TAIL, "Genomic position CSS exponent for tail", GEN_POS_CSS_EXPONENT_TAIL);
        configBuilder.addDecimal(CSS_EXPONENT_GEN_POS_CUTOFF, "Genomic position CSS exponent cutoff", GEN_POS_CSS_EXPONENT_CUTOFF);
        configBuilder.addInteger(GEN_POS_BUCKET_SIZE_CFG, GEN_POS_BUCKET_SIZE_DESC, GEN_POS_BUCKET_SIZE);
        configBuilder.addInteger(GEN_POS_MAX_SAMPLE_COUNT_CFG, GEN_POS_MAX_SAMPLE_COUNT_DESC, GEN_POS_MAX_SAMPLE_COUNT);
        configBuilder.addFlag(WRITE_SNV_SIMILARITIES, "Write SNV-96 CSS to file");
        configBuilder.addFlag(WRITE_GEN_POS_SIMILARITIES, "Write genomic position CSS to file");
        configBuilder.addFlag(WRITE_GEN_POS_CSS, "Write gen-pos CSS to file");
    }

    public CategoryType categoryType() { return SNV; }

    public void close()
    {
        closeBufferedWriter(mGenPosCohortCssWriter);
    }

    @Override
    public boolean loadData()
    {
        // first load data
        if(mConfig.RefSnvCountsFile.isEmpty() && mConfig.RefSigContributionFile.isEmpty() && mConfig.RefSnvCancerPosFreqFile.isEmpty())
            return false;

        if(!loadRefSignaturePercentileData(
                mConfig.RefSigContributionFile, mSigContributions.getRefCancerSigContribPercentiles(), mRefCancerSnvCountPercentiles))
        {
            return false;
        }

        mRefSampleSnv96Counts = loadSampleCountsFromFile(mConfig.RefSnvCountsFile, mRefSampleSnv96CountsIndex);

        mRefCancerGenPosCounts = loadRefSampleCounts(mConfig.RefSnvCancerPosFreqFile, mRefGenPosCancerTypes, Lists.newArrayList());
        mRefSampleGenPosCounts = loadSampleMatrixData(mConfig.RefSnvSamplePosFreqFile, mRefSampleGenPosCountsIndex);

        if(mRefSampleSnv96Counts == null || mRefCancerGenPosCounts == null)
        {
            CUP_LOGGER.error("invalid somatic matrix data: SNV-96({}) GenPos({})",
                    mRefSampleSnv96Counts != null, mRefCancerGenPosCounts != null);
            return false;
        }

        if(!mConfig.TestRefData)
        {
            if(!loadSampleCounts())
                return false;
        }
        else
        {
            mSampleSnv96Counts = mRefSampleSnv96Counts;
            mSampleSnv96CountsIndex.putAll(mRefSampleSnv96CountsIndex);
            mSampleGenPosCounts = mRefSampleGenPosCounts;
            mSampleGenPosCountsIndex.putAll(mRefSampleGenPosCountsIndex);
        }

        if(!mSigContributions.loadSigContributions(mSampleSnv96Counts, mSampleSnv96CountsIndex))
            return false;

        // record SNV totals, prior to any noise adjustments
        for(Map.Entry<String,Integer> entry : mRefSampleSnv96CountsIndex.entrySet())
        {
            String sampleId = entry.getKey();
            final double[] sampleCounts = mRefSampleSnv96Counts.getRow(entry.getValue());
            int snvTotal = (int) sumVector(sampleCounts);
            mSampleSnvTotals.put(sampleId, snvTotal);
        }

        if(!mConfig.TestRefData)
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

            if(!mConfig.TestRefData)
                NoiseRefCache.applyNoise(mSampleSnv96Counts, noiseAdjustments, noiseAllocation);
        }

        return true;
    }

    private boolean loadSampleCounts()
    {
        int sampleCount = mSampleDataCache.SampleDataList.size();

        final Map<String,Integer> triNucBucketNameMap = Maps.newHashMap();
        populateBucketMap(triNucBucketNameMap);

        mSampleSnv96Counts = new Matrix(sampleCount, triNucBucketNameMap.size());
        mSampleGenPosCounts = new Matrix(sampleCount, mPosFrequencies.getBucketCount());

        for(int i = 0; i < sampleCount; ++i)
        {
            final String sampleId = mSampleDataCache.SampleIds.get(i);

            List<SomaticVariant> somaticVariants = Lists.newArrayList();

            // load from VCF or database
            if(mConfig.DbAccess != null)
            {
                somaticVariants.addAll(loadSomaticVariants(sampleId, mConfig.DbAccess));
            }
            else
            {
                String purpleDir = mConfig.getPurpleDataDir(sampleId);
                final String somaticVcfFile = purpleSomaticVcfFile(purpleDir, sampleId);
                somaticVariants.addAll(SomaticDataLoader.loadSomaticVariantsFromVcf(somaticVcfFile, Lists.newArrayList(SNP)));
            }

            final double[] triNucCounts = extractTrinucleotideCounts(somaticVariants, triNucBucketNameMap);

            mSampleSnv96CountsIndex.put(sampleId, i);
            mSampleSnv96Counts.setRow(i, triNucCounts);

            // mSampleSnv96Counts = convertSomaticVariantsToSnvCounts(sampleId, somaticVariants, mSampleSnv96CountsIndex);

            AidApobecStatus aidApobecStatus = AidApobecStatus.FALSE_ONLY;

            mPosFrequencies.clear();
            extractPositionFrequencyCounts(somaticVariants, mPosFrequencies, aidApobecStatus);

            mSampleGenPosCounts.setRow(i, mPosFrequencies.getCounts());
            mSampleGenPosCountsIndex.put(sampleId, i);

            /* don't bother using these sig counts files, just regenerate on the fly
            final String snvCountsFile = !mConfig.SampleSnvCountsFile.isEmpty() ?
                    mConfig.SampleSnvCountsFile : mConfig.SampleDataDir + sampleId + ".sig.snv_counts.csv";

            final String snvPosFreqFile = !mConfig.SampleSnvPosFreqFile.isEmpty() ?
                    mConfig.SampleSnvPosFreqFile : mConfig.SampleDataDir + sampleId + ".sig.pos_freq_counts.csv";

            mSampleSnv96Counts = loadSampleCountsFromFile(snvCountsFile, mSampleSnv96CountsIndex);
            mSampleGenPosCounts = loadSampleMatrixData(snvPosFreqFile, mSampleGenPosCountsIndex);
            */
        }

        // no support for cohort matrix files for non-ref cohorts

        /*
        if(!mConfig.SampleSnvCountsFile.isEmpty() && !mConfig.SampleSnvPosFreqFile.isEmpty())
        {
            mSampleSnv96Counts = loadSampleCountsFromFile(mConfig.SampleSnvCountsFile, mSampleSnv96CountsIndex);

            if(mSampleSnv96Counts == null)
            {
                CUP_LOGGER.error("missing file: {}", mConfig.SampleSnvCountsFile);
            }

            mSampleGenPosCounts = loadSampleMatrixData(mConfig.SampleSnvPosFreqFile, mSampleGenPosCountsIndex);

            if(mSampleSnv96Counts == null)
            {
                CUP_LOGGER.error("missing file: {}", mConfig.SampleSnvPosFreqFile);
            }

            return mSampleSnv96Counts != null && mSampleGenPosCounts != null;
        }
        */

        return true;
    }

    @Override
    public boolean processSample(final SampleData sample, final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        if(mRefSampleSnv96Counts == null)
            return false;

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

        return true;
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

        addCohortGenPosCssResults(sample, sampleCounts, sampleTotal, results, similarities);
    }

    private void addCohortGenPosCssResults(
            final SampleData sample, final double[] sampleCounts, double snvTotal,
            final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        // first run CSS against cancer cohorts
        int refCancerCount = mRefCancerGenPosCounts.Rows;
        double maxCssScore = 0;

        final Map<String,Double> cancerCssTotals = Maps.newHashMap();

        int maxSampleCount = mPosFrequencies.getMaxSampleCount();

        for(int i = 0; i < refCancerCount; ++i)
        {
            final String refCancerType = mRefGenPosCancerTypes.get(i);

            if(!isKnownCancerType(refCancerType))
                continue;

            if(!checkIsValidCancerType(sample, refCancerType, cancerCssTotals))
                continue;

            boolean matchesCancerType = sample.cancerType().equals(refCancerType);

            double adjustMultiplier = snvTotal > maxSampleCount ? maxSampleCount / snvTotal : 1;

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

        results.add(new SampleResult(
                sample.Id, SNV, CLASSIFIER, GENOMIC_POSITION_COHORT.toString(), String.format("%.4g", totalWeightedCss), cancerCssTotals));
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
