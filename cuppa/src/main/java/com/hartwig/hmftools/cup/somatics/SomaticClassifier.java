package com.hartwig.hmftools.cup.somatics;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;
import static java.lang.Math.sqrt;

import static com.hartwig.hmftools.common.stats.CosineSimilarity.calcCosineSim;
import static com.hartwig.hmftools.common.stats.Percentiles.getPercentile;
import static com.hartwig.hmftools.common.sigs.VectorUtils.sumVector;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CategoryType.CLASSIFIER;
import static com.hartwig.hmftools.cup.common.CategoryType.SAMPLE_TRAIT;
import static com.hartwig.hmftools.cup.common.CategoryType.SNV;
import static com.hartwig.hmftools.cup.common.ClassifierType.SNV_96_PAIRWISE_SIMILARITY;
import static com.hartwig.hmftools.cup.common.ClassifierType.GENOMIC_POSITION_SIMILARITY;
import static com.hartwig.hmftools.cup.common.CupCalcs.adjustRefCounts;
import static com.hartwig.hmftools.cup.common.CupCalcs.calcPercentilePrevalence;
import static com.hartwig.hmftools.cup.common.CupCalcs.convertToPercentages;
import static com.hartwig.hmftools.cup.common.CupCalcs.fillMissingCancerTypeValues;
import static com.hartwig.hmftools.cup.common.CupConstants.CSS_SIMILARITY_CUTOFF;
import static com.hartwig.hmftools.cup.common.CupConstants.CSS_SIMILARITY_MAX_MATCHES;
import static com.hartwig.hmftools.cup.common.CupConstants.SNV_CSS_DIFF_EXPONENT;
import static com.hartwig.hmftools.cup.common.CupConstants.SNV_CSS_THRESHOLD;
import static com.hartwig.hmftools.cup.common.CupConstants.SNV_POS_FREQ_CSS_THRESHOLD;
import static com.hartwig.hmftools.cup.common.CupConstants.SNV_POS_FREQ_DIFF_EXPONENT;
import static com.hartwig.hmftools.cup.common.ResultType.LIKELIHOOD;
import static com.hartwig.hmftools.cup.common.ResultType.PERCENTILE;
import static com.hartwig.hmftools.cup.common.SampleData.isKnownCancerType;
import static com.hartwig.hmftools.cup.common.SampleResult.checkIsValidCancerType;
import static com.hartwig.hmftools.cup.common.SampleSimilarity.recordCssSimilarity;
import static com.hartwig.hmftools.cup.somatics.RefSomatics.convertSignatureName;
import static com.hartwig.hmftools.cup.somatics.RefSomatics.populateReportableSignatures;
import static com.hartwig.hmftools.cup.somatics.SomaticDataLoader.loadRefSampleCounts;
import static com.hartwig.hmftools.cup.somatics.SomaticDataLoader.loadRefSignaturePercentileData;
import static com.hartwig.hmftools.cup.somatics.SomaticDataLoader.loadSampleCountsFromFile;
import static com.hartwig.hmftools.cup.somatics.SomaticDataLoader.loadSamplePosFreqFromFile;
import static com.hartwig.hmftools.cup.somatics.SomaticDataLoader.loadSigContribsFromCohortFile;
import static com.hartwig.hmftools.cup.somatics.SomaticDataLoader.loadSigContribsFromDatabase;

import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.common.sigs.SignatureAllocation;
import com.hartwig.hmftools.common.sigs.SignatureAllocationFile;
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
    private final Map<String,Map<String,double[]>> mRefCancerSigContribPercentiles;
    private final Map<String,double[]> mRefCancerSnvCountPercentiles;

    private Matrix mRefCancerSnvPosFrequencies;
    private final List<String> mRefSnvPosFreqCancerTypes;

    private Matrix mRefSamplePosFrequencies;
    private final Map<String,Integer> mRefSamplePosFreqIndex;

    private Matrix mSampleCounts;
    private final Map<String,Integer> mSampleCountsIndex;

    private final Map<String,Map<String,Double>> mSampleSigContributions;

    private Matrix mSamplePosFrequencies;
    private final Map<String,Integer> mSamplePosFreqIndex;

    private boolean mIsValid;

    private final double mMaxCssAdjustFactor;
    private static final int SNV_POS_FREQ_SNV_TOTAL_THRESHOLD = 20000;

    //private static final double MAX_CSS_ADJUST_FACTOR = 4;
    public static final String MAX_CSS_ADJUST_FACTOR = "css_max_factor";

    public SomaticClassifier(final CuppaConfig config, final SampleDataCache sampleDataCache, final CommandLine cmd)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mSampleCounts = null;
        mRefCancerSnvPosFrequencies = null;
        mRefSamplePosFrequencies = null;
        mSampleSigContributions = Maps.newHashMap();
        mSampleCountsIndex = Maps.newHashMap();
        mSamplePosFreqIndex = Maps.newHashMap();

        mRefSampleCounts = null;
        mRefSampleNames = Lists.newArrayList();
        mRefCancerSigContribPercentiles = Maps.newHashMap();
        mRefCancerSnvCountPercentiles = Maps.newHashMap();
        mRefSnvPosFreqCancerTypes = Lists.newArrayList();
        mRefSamplePosFreqIndex = Maps.newHashMap();

        mMaxCssAdjustFactor = cmd != null ? Double.parseDouble(cmd.getOptionValue(MAX_CSS_ADJUST_FACTOR, "0")) : 0;

        mIsValid = true;

        if(mConfig.RefSnvCountsFile.isEmpty() && mConfig.RefSigContributionFile.isEmpty() && mConfig.RefSnvCancerPosFreqFile.isEmpty())
            return;

        populateReportableSignatures();

        loadRefSignaturePercentileData(mConfig.RefSigContributionFile, mRefCancerSigContribPercentiles, mRefCancerSnvCountPercentiles);
        mRefSampleCounts = loadRefSampleCounts(mConfig.RefSnvCountsFile, mRefSampleNames, Lists.newArrayList("BucketName"));

        mRefCancerSnvPosFrequencies = loadRefSampleCounts(mConfig.RefSnvCancerPosFreqFile, mRefSnvPosFreqCancerTypes, Lists.newArrayList());
        mRefSamplePosFrequencies = loadSamplePosFreqFromFile(mConfig.RefSnvSamplePosFreqFile, mRefSamplePosFreqIndex);

        if(mRefSampleCounts == null || mRefCancerSnvPosFrequencies == null)
            mIsValid = false;

        mIsValid &= loadSampleCounts();
        mIsValid &= loadSigContributions();
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(MAX_CSS_ADJUST_FACTOR, true, "Max CSS adustment factor");
    }

    public CategoryType categoryType() { return SNV; }
    public boolean isValid() { return mIsValid; }
    public void close() {}

    private boolean loadSampleCounts()
    {
        if(!mConfig.SampleSnvCountsFile.isEmpty() && !mConfig.SampleSnvPosFreqFile.isEmpty())
        {
            if(mConfig.SampleSnvCountsFile.equals(mConfig.RefSnvCountsFile))
            {
                mSampleCounts = mRefSampleCounts;

                for(int i = 0; i < mRefSampleNames.size(); ++i)
                {
                    mSampleCountsIndex.put(mRefSampleNames.get(i), i);
                }
            }
            else
            {
                mSampleCounts = loadSampleCountsFromFile(mConfig.SampleSnvCountsFile, mSampleCountsIndex);
            }

            if(mConfig.SampleSnvPosFreqFile.equals(mConfig.RefSnvSamplePosFreqFile))
            {
                mSamplePosFrequencies = mRefSamplePosFrequencies;
                mSamplePosFreqIndex.putAll(mRefSamplePosFreqIndex);
            }
            else
            {
                mSamplePosFrequencies = loadSamplePosFreqFromFile(mConfig.SampleSnvPosFreqFile, mSamplePosFreqIndex);
            }

            return mSampleCounts != null && mSamplePosFrequencies != null;
        }

        if(mSampleDataCache.isSingleSample())
        {
            final String sampleId = mSampleDataCache.SampleIds.get(0);

            /*
            if(!mConfig.SampleSomaticVcf.isEmpty())
            {
                final List<SomaticVariant> somaticVariants = loadSomaticVariants(sampleId, mConfig.SampleSomaticVcf);
                populateSomaticCounts(somaticVariants);
            }
            */

            final String snvCountsFile = !mConfig.SampleSnvCountsFile.isEmpty() ?
                    mConfig.SampleSnvCountsFile : mConfig.SampleDataDir + sampleId + ".sig.snv_counts.csv";

            final String snvPosFreqFile = !mConfig.SampleSnvPosFreqFile.isEmpty() ?
                    mConfig.SampleSnvPosFreqFile : mConfig.SampleDataDir + sampleId + ".sig.pos_freq_counts.csv";

            mSampleCounts = loadSampleCountsFromFile(snvCountsFile, mSampleCountsIndex);
            mSamplePosFrequencies = loadSamplePosFreqFromFile(snvPosFreqFile, mSamplePosFreqIndex);

            return mSampleCounts != null && mSamplePosFrequencies != null;
        }
        else if(mConfig.DbAccess != null)
        {
            CUP_LOGGER.error("somatic variants from DB not supported");

            /*
            final String sampleId = mSampleDataCache.SampleIds.get(0);
            final List<SomaticVariant> somaticVariants = loadSomaticVariants(sampleId, mConfig.DbAccess);
            populateSomaticCounts(somaticVariants);

            return mSampleCounts != null && mSamplePosFrequencies != null;
            */

            return false;
        }
        else
        {
            CUP_LOGGER.error("no sample SNV count source specified");
            return false;
        }
    }

    private boolean loadSigContributions()
    {
        if(!mConfig.SampleSigContribFile.isEmpty())
        {
            CUP_LOGGER.info("loading SNV sig contributions from file({})", mConfig.SampleSigContribFile);
            return loadSigContribsFromCohortFile(mConfig.SampleSigContribFile, mSampleSigContributions);
        }
        else if(mConfig.DbAccess != null)
        {
            CUP_LOGGER.info("loading SNV sig contributions from database");
            return loadSigContribsFromDatabase(mConfig.DbAccess, mSampleDataCache.SampleIds, mSampleSigContributions);
        }

        final String sampleId = mSampleDataCache.SampleIds.get(0);
        final String sigAllocFile = SignatureAllocationFile.generateFilename(mConfig.SampleDataDir, sampleId);

        try
        {
            final List<SignatureAllocation> sigAllocations = SignatureAllocationFile.read(sigAllocFile);
            Map<String,Double> sigContribs = Maps.newHashMap();
            for(final SignatureAllocation sigAllocation : sigAllocations)
            {
                final String sigName = convertSignatureName(sigAllocation.signature());
                sigContribs.put(sigName, sigAllocation.allocation());
            }

            mSampleSigContributions.put(sampleId, sigContribs);
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("sample({}) failed to load sig allocations file({}): {}",
                    sampleId, sigAllocFile, e.toString());
            return false;
        }

        return true;
    }

    public void processSample(final SampleData sample, final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        if(!mIsValid || mRefSampleCounts == null)
            return;

        Integer sampleCountsIndex = mSampleCountsIndex.get(sample.Id);

        if(sampleCountsIndex == null)
        {
            CUP_LOGGER.info("sample({}) has no SNV data", sample.Id);
            return;
        }

        final double[] sampleCounts = mSampleCounts.getCol(sampleCountsIndex);
        int snvTotal = (int)sumVector(sampleCounts);

        addCssResults(sample, sampleCounts, snvTotal, results, similarities);
        addPosFreqCssResults(sample, results, similarities);

        addSigContributionResults(sample, results);

        // add a percentile result
        final Map<String, Double> cancerTypeValues = Maps.newHashMap();

        for(Map.Entry<String, double[]> cancerPercentiles : mRefCancerSnvCountPercentiles.entrySet())
        {
            final String cancerType = cancerPercentiles.getKey();

            if(!isKnownCancerType(cancerType))
                continue;

            double percentile = getPercentile(cancerPercentiles.getValue(), snvTotal, true);
            cancerTypeValues.put(cancerType, percentile);
        }

        SampleResult result = new SampleResult(
                sample.Id, SAMPLE_TRAIT, PERCENTILE, "SNV_COUNT", snvTotal, cancerTypeValues);

        results.add(result);

        int cancerTypeCount = mSampleDataCache.RefCancerSampleData.size();
        int cancerSampleCount = sample.isRefSample() ? mSampleDataCache.getCancerSampleCount(sample.cancerType()) : 0;

        final Map<String,Double> cancerPrevsLow = calcPercentilePrevalence(
                sample, cancerSampleCount, cancerTypeCount, mRefCancerSnvCountPercentiles, snvTotal,  true);

        results.add(new SampleResult(sample.Id, SNV, LIKELIHOOD, "SNV_COUNT_LOW", snvTotal, cancerPrevsLow));

        final Map<String,Double> cancerPrevsHigh = calcPercentilePrevalence(
                sample, cancerSampleCount, cancerTypeCount, mRefCancerSnvCountPercentiles, snvTotal, false);

        results.add(new SampleResult(sample.Id, SNV, LIKELIHOOD, "SNV_COUNT_HIGH", snvTotal, cancerPrevsHigh));
    }

    private void addCssResults(
            final SampleData sample, final double[] sampleCounts, int snvTotal,
            final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        int refSampleCount = mRefSampleCounts.Cols;

        final List<SampleSimilarity> topMatches = Lists.newArrayList();

        final Map<String,Double> cancerCssTotals = Maps.newHashMap();

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

            if(mConfig.WriteSimilarities)
            {
                recordCssSimilarity(
                        topMatches, sample.Id, refSampleId, css, SNV_96_PAIRWISE_SIMILARITY.toString(),
                        CSS_SIMILARITY_MAX_MATCHES, CSS_SIMILARITY_CUTOFF);
            }

            if(!isKnownCancerType(refCancerType))
                continue;

            maxCssScore = max(css, maxCssScore);

            double cssWeight = pow(SNV_CSS_DIFF_EXPONENT, -100 * (1 - css));

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

        if(mMaxCssAdjustFactor > 0 && totalCss > 0)
        {
            double adjustedFactor = mMaxCssAdjustFactor > 0 ? pow(maxCssScore, mMaxCssAdjustFactor) : 0;

            for(Map.Entry<String,Double> entry : cancerCssTotals.entrySet())
            {
                double adjCancerScore = pow(entry.getValue(), adjustedFactor);
                cancerCssTotals.put(entry.getKey(), adjCancerScore);
            }

            convertToPercentages(cancerCssTotals);
        }

        results.add(new SampleResult(
                sample.Id, CLASSIFIER, LIKELIHOOD, SNV_96_PAIRWISE_SIMILARITY.toString(), String.format("%.4g", totalCss), cancerCssTotals));

        // for non-ref cohorts, also report closest matches from amongst these
        if(mConfig.WriteSimilarities && mSampleDataCache.isMultiSampleNonRef())
        {
            for(Map.Entry<String,Integer> entry : mSampleCountsIndex.entrySet())
            {
                final String nonRefSampleId = entry.getKey();

                if(nonRefSampleId.equals(sample.Id))
                    continue;

                final double[] otherSampleCounts = mSampleCounts.getCol(entry.getValue());

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

        final double[] sampleCounts = mSamplePosFrequencies.getCol(sampleCountsIndex);
        double sampleTotal = sumVector(sampleCounts);

        // first run CSS against cancer cohorts
        int refCancerCount = mRefCancerSnvPosFrequencies.Cols;

        final Map<String,Double> cancerCssTotals = Maps.newHashMap();

        for(int i = 0; i < refCancerCount; ++i)
        {
            final String refCancerType = mRefSnvPosFreqCancerTypes.get(i);

            if(!isKnownCancerType(refCancerType))
                continue;

            if(!checkIsValidCancerType(sample, refCancerType, cancerCssTotals))
                continue;

            boolean matchesCancerType = sample.cancerType().equals(refCancerType);

            double adjustMultiplier = sampleTotal > SNV_POS_FREQ_SNV_TOTAL_THRESHOLD ? SNV_POS_FREQ_SNV_TOTAL_THRESHOLD / sampleTotal : 1;

            final double[] refPosFreqs = sample.isRefSample() && matchesCancerType ?
                    adjustRefCounts(mRefCancerSnvPosFrequencies.getCol(i), sampleCounts, adjustMultiplier) : mRefCancerSnvPosFrequencies.getCol(i);

            double css = calcCosineSim(sampleCounts, refPosFreqs);

            if(css < SNV_POS_FREQ_CSS_THRESHOLD)
                continue;

            double cssWeight = pow(SNV_POS_FREQ_DIFF_EXPONENT, -100 * (1 - css));

            double weightedCss = css * cssWeight;

            Double total = cancerCssTotals.get(refCancerType);

            if(total == null)
                cancerCssTotals.put(refCancerType, weightedCss);
            else
                cancerCssTotals.put(refCancerType, total + weightedCss);
        }

        double totalCss = cancerCssTotals.values().stream().mapToDouble(x -> x).sum();

        convertToPercentages(cancerCssTotals);

        results.add(new SampleResult(
                sample.Id, CLASSIFIER, LIKELIHOOD, GENOMIC_POSITION_SIMILARITY.toString(), String.format("%.4g", totalCss), cancerCssTotals));

        // then run pairwise analysis if similarities are being analysed
        // addSnvPosSimilarities(sample, sampleCounts, similarities);
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

            if(mConfig.WriteSimilarities)
            {
                recordCssSimilarity(
                        topMatches, sample.Id, refSampleId, css, GENOMIC_POSITION_SIMILARITY.toString(),
                        CSS_SIMILARITY_MAX_MATCHES, CSS_SIMILARITY_CUTOFF);
            }
        }

        similarities.addAll(topMatches);
    }

    private static final String SIG_NAME_2 = "Sig2";
    private static final String SIG_NAME_13 = "Sig13";

    private static final String signatureDisplayName(final String sigName)
    {
        final String displayName = RefSomatics.REPORTABLE_SIGS.get(sigName);
        return displayName != null ? displayName : "UNKNOWN";
    }

    private void addSigContributionResults(final SampleData sample, final List<SampleResult> results)
    {
        final Map<String,Double> sampleSigContribs = mSampleSigContributions.get(sample.Id);

        if(sampleSigContribs == null)
        {
            CUP_LOGGER.debug("sample({}) sig contributions not found", sample.Id);
            return;
        }

        // report on every one of the designated set

        for(final String sigName : RefSomatics.REPORTABLE_SIGS.keySet())
        {
            double sampleSigContrib = sampleSigContribs.containsKey(sigName) ? sampleSigContribs.get(sigName) : 0;

            // report the AID/APOBEC sigs 2 & 13 together
            if(sigName.equalsIgnoreCase(SIG_NAME_2))
            {
                Double otherAlloc = sampleSigContribs.get(SIG_NAME_13);
                if(otherAlloc != null)
                    sampleSigContrib += otherAlloc;
            }
            else if(sigName.equalsIgnoreCase(SIG_NAME_13))
            {
                continue;
            }

            Map<String, Double> cancerResults = Maps.newHashMap();

            for(Map.Entry<String,Map<String,double[]>> cancerContribs : mRefCancerSigContribPercentiles.entrySet())
            {
                final String cancerType = cancerContribs.getKey();
                final double[] refSigPercentiles = cancerContribs.getValue().get(sigName);

                if(refSigPercentiles == null)
                {
                    // CUP_LOGGER.debug("missing sig({}) data for cancerType({})", sigName, cancerType);
                    cancerResults.put(cancerType, 0.0);
                }
                else
                {
                    double percentile = getPercentile(refSigPercentiles, sampleSigContrib, true);
                    cancerResults.put(cancerType, percentile);
                }
            }

            results.add(new SampleResult(sample.Id, SNV, PERCENTILE, signatureDisplayName(sigName), round(sampleSigContrib), cancerResults));
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
        mSampleCounts = new Matrix(snvCounts.get(0).length, snvCounts.size());

        for(int i = 0; i < snvCounts.size(); ++i)
        {
            mSampleCounts.setCol(i, snvCounts.get(i));
            mSampleCountsIndex.put(sampleIds.get(i), i);
        }

        mSamplePosFrequencies = new Matrix(posFreqCounts.get(0).length, posFreqCounts.size());

        for(int i = 0; i < posFreqCounts.size(); ++i)
        {
            mSamplePosFrequencies.setCol(i, posFreqCounts.get(i));
            mSamplePosFreqIndex.put(sampleIds.get(i), i);
        }
    }

}
