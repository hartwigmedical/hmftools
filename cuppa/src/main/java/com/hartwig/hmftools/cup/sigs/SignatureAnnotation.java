package com.hartwig.hmftools.cup.sigs;

import static java.lang.Math.max;
import static java.lang.Math.pow;
import static java.lang.Math.round;
import static java.lang.Math.sqrt;

import static com.hartwig.hmftools.common.sigs.CosineSimilarity.calcCosineSim;
import static com.hartwig.hmftools.common.sigs.Percentiles.getPercentile;
import static com.hartwig.hmftools.common.sigs.VectorUtils.copyVector;
import static com.hartwig.hmftools.common.sigs.VectorUtils.sumVector;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CategoryType.CLASSIFIER;
import static com.hartwig.hmftools.cup.common.CategoryType.SAMPLE_TRAIT;
import static com.hartwig.hmftools.cup.common.CategoryType.SNV_SIG;
import static com.hartwig.hmftools.cup.common.ClassifierType.SNV_COUNT_CSS;
import static com.hartwig.hmftools.cup.common.ClassifierType.SNV_POS_FREQ_CSS;
import static com.hartwig.hmftools.cup.common.CupConstants.SNV_CSS_THRESHOLD;
import static com.hartwig.hmftools.cup.common.CupConstants.SNV_POS_FREQ_CSS_THRESHOLD;
import static com.hartwig.hmftools.cup.sigs.SignatureDataLoader.loadRefSampleCounts;
import static com.hartwig.hmftools.cup.sigs.SignatureDataLoader.loadRefSigContribPercentiles;
import static com.hartwig.hmftools.cup.sigs.SignatureDataLoader.loadRefSnvPosFrequences;
import static com.hartwig.hmftools.cup.sigs.SignatureDataLoader.loadSampleCountsFromCohortFile;
import static com.hartwig.hmftools.cup.sigs.SignatureDataLoader.loadSamplePosFreqFromCohortFile;
import static com.hartwig.hmftools.cup.sigs.SignatureDataLoader.loadSigContribsFromCohortFile;
import static com.hartwig.hmftools.cup.sigs.SignatureDataLoader.loadSigContribsFromDatabase;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sigs.SigMatrix;
import com.hartwig.hmftools.cup.SampleAnalyserConfig;
import com.hartwig.hmftools.cup.common.ClassifierType;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;

public class SignatureAnnotation
{
    private final SampleAnalyserConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private SigMatrix mRefSampleCounts;
    private final List<String> mRefSampleNames;
    private final Map<String,Map<String,double[]>> mRefCancerSigContribPercentiles;

    private SigMatrix mRefSnvPosFrequencies;
    private final List<String> mRefSnvPosFreqCancerTypes;

    private SigMatrix mSampleCounts;
    private final Map<String,Integer> mSampleCountsIndex;

    private final Map<String,Map<String,Double>> mSampleSigContributions;

    private SigMatrix mSamplePosFrequencies;
    private final Map<String,Integer> mSamplePosFreqIndex;

    private static final int SNV_POS_FREQ_SNV_TOTAL_THRESHOLD = 20000;

    public SignatureAnnotation(final SampleAnalyserConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mSampleCounts = null;
        mRefSnvPosFrequencies = null;
        mSampleSigContributions = Maps.newHashMap();
        mSampleCountsIndex = Maps.newHashMap();
        mSamplePosFreqIndex = Maps.newHashMap();

        mRefSampleCounts = null;
        mRefSampleNames = Lists.newArrayList();
        mRefCancerSigContribPercentiles = Maps.newHashMap();

        loadRefSigContribPercentiles(mConfig.RefSigContribData, mRefCancerSigContribPercentiles);
        mRefSampleCounts = loadRefSampleCounts(mConfig.RefSnvCountsFile, mRefSampleNames);

        mRefSnvPosFreqCancerTypes = Lists.newArrayList();
        mRefSnvPosFrequencies = loadRefSnvPosFrequences(mConfig.RefSnvPosFreqFile, mRefSnvPosFreqCancerTypes);

        loadSampleCounts();
        loadSigContributions();
    }

    private boolean loadSampleCounts()
    {
        if(!mConfig.SampleSnvCountsFile.isEmpty())
        {
            mSampleCounts = loadSampleCountsFromCohortFile(mConfig.SampleSnvCountsFile, mSampleCountsIndex);
            mSamplePosFrequencies = loadSamplePosFreqFromCohortFile(mConfig.SampleSnvPosFreqFile, mSamplePosFreqIndex);
        }
        else if(mConfig.DbAccess != null)
        {
            CUP_LOGGER.error("retrieval of sample SNV counts unsupported at present");
            return false;
        }

        return true;
    }

    private void loadSigContributions()
    {
        if(!mConfig.SampleSigContribFile.isEmpty())
        {
            loadSigContribsFromCohortFile(mConfig.SampleSigContribFile, mSampleSigContributions);
        }
        else if(mConfig.DbAccess != null)
        {
            loadSigContribsFromDatabase(mConfig.DbAccess, mSampleDataCache.SampleIds, mSampleSigContributions);
        }
    }

    public int getSampleSnvCount(final String sampleId)
    {
        Integer sampleCountsIndex = mSampleCountsIndex.get(sampleId);

        if(sampleCountsIndex == null)
        {
            CUP_LOGGER.info("sample({}) has no SNV data", sampleId);
            return 0;
        }

        final double[] sampleCounts = mSampleCounts.getCol(sampleCountsIndex);
        return (int)sumVector(sampleCounts);
    }

    public List<SampleResult> processSample(final SampleData sample)
    {
        final List<SampleResult> results = Lists.newArrayList();

        Integer sampleCountsIndex = mSampleCountsIndex.get(sample.Id);

        if(sampleCountsIndex == null)
        {
            CUP_LOGGER.info("sample({}) has no SNV data", sample.Id);
            return results;
        }

        final double[] sampleCounts = mSampleCounts.getCol(sampleCountsIndex);
        int snvTotal = (int)sumVector(sampleCounts);

        if(mConfig.runCategory(CLASSIFIER))
        {
            addCssResults(sample, sampleCounts, results);
            addPosFreqCssResults(sample, results);
        }

        if(mConfig.runCategory(SNV_SIG))
        {
            addSigContributionResults(sample, snvTotal, results);
        }

        return results;
    }

    private void addCssResults(final SampleData sample, final double[] sampleCounts, final List<SampleResult> results)
    {
        int refSampleCount = mRefSampleCounts.Cols;

        final Map<String,Double> cancerCssTotals = Maps.newHashMap();

        for(int s = 0; s < refSampleCount; ++s)
        {
            final String refSampleId = mRefSampleNames.get(s);

            if(refSampleId.equals(sample.Id))
                continue;

            final double[] otherSampleCounts = mRefSampleCounts.getCol(s);

            double css = calcCosineSim(sampleCounts, otherSampleCounts);

            if(css < SNV_CSS_THRESHOLD)
                continue;

            final String refCancerType = mSampleDataCache.RefSampleCancerTypeMap.get(refSampleId);

            double cssWeight = pow(2, -100 * (1 - css));

            int cancerTypeCount = mSampleDataCache.RefCancerSampleData.get(refCancerType).size();
            double weightedCss = css * cssWeight / sqrt(cancerTypeCount);

            Double total = cancerCssTotals.get(refCancerType);

            if(total == null)
                cancerCssTotals.put(refCancerType, weightedCss);
            else
                cancerCssTotals.put(refCancerType, total + weightedCss);
        }

        double totalCss = cancerCssTotals.values().stream().mapToDouble(x -> x).sum();

        for(Map.Entry<String,Double> entry : cancerCssTotals.entrySet())
        {
            cancerCssTotals.put(entry.getKey(), entry.getValue() / totalCss);
        }

        results.add(new SampleResult(sample.Id, CLASSIFIER, ClassifierType.displayString(SNV_COUNT_CSS), String.format("%.6f", totalCss), cancerCssTotals));
    }

    private void addPosFreqCssResults(final SampleData sample, final List<SampleResult> results)
    {
        Integer sampleCountsIndex = mSamplePosFreqIndex.get(sample.Id);

        if(sampleCountsIndex == null)
        {
            CUP_LOGGER.info("sample({}) has no SNV pos-freq data", sample.Id);
            return;
        }

        final double[] sampleCounts = mSamplePosFrequencies.getCol(sampleCountsIndex);
        double sampleTotal = sumVector(sampleCounts);

        int refCancerCount = mRefSnvPosFrequencies.Cols;

        final Map<String,Double> cancerCssTotals = Maps.newHashMap();

        for(int i = 0; i < refCancerCount; ++i)
        {
            final String refCancerType = mRefSnvPosFreqCancerTypes.get(i);
            boolean matchesCancerType = sample.CancerType.equals(refCancerType);

            final double[] refPosFreqs = mConfig.AdjustSnvPosFreqCounts && matchesCancerType ?
                    adjustRefPosFreqCounts(mRefSnvPosFrequencies.getCol(i), sampleCounts, sampleTotal) : mRefSnvPosFrequencies.getCol(i);

            double css = calcCosineSim(sampleCounts, refPosFreqs);

            if(css < SNV_POS_FREQ_CSS_THRESHOLD)
                continue;

            double cssWeight = pow(10, -100 * (1 - css));

            double weightedCss = css * cssWeight;

            Double total = cancerCssTotals.get(refCancerType);

            if(total == null)
                cancerCssTotals.put(refCancerType, weightedCss);
            else
                cancerCssTotals.put(refCancerType, total + weightedCss);
        }

        double totalCss = cancerCssTotals.values().stream().mapToDouble(x -> x).sum();

        for(Map.Entry<String,Double> entry : cancerCssTotals.entrySet())
        {
            cancerCssTotals.put(entry.getKey(), entry.getValue() / totalCss);
        }

        results.add(new SampleResult(sample.Id, CLASSIFIER, ClassifierType.displayString(SNV_POS_FREQ_CSS), String.format("%.6f", totalCss), cancerCssTotals));
    }

    private double[] adjustRefPosFreqCounts(final double[] refPosFreqs, final double[] sampleCounts, final double sampleTotal)
    {
        double adjustMultiplier = sampleTotal > SNV_POS_FREQ_SNV_TOTAL_THRESHOLD ? SNV_POS_FREQ_SNV_TOTAL_THRESHOLD / sampleTotal : 1;

        double[] adjustedCounts = new double[refPosFreqs.length];

        for(int b = 0; b < refPosFreqs.length; ++b)
        {
            adjustedCounts[b] = max(refPosFreqs[b] - (sampleCounts[b] * adjustMultiplier), 0);
        }

        return adjustedCounts;
    }

    private static final List<String> REPORTABLE_SIGS =
            Lists.newArrayList("Sig1", "Sig2", "Sig4", "Sig6", "Sig7", "Sig10", "Sig11", "Sig13", "Sig17");

    private static final String SIG_NAME_2 = "Sig2";
    private static final String SIG_NAME_13 = "Sig13";

    private void addSigContributionResults(final SampleData sample, int snvTotal, final List<SampleResult> results)
    {
        final Map<String,Double> sampleSigContribs = mSampleSigContributions.get(sample.Id);

        if(sampleSigContribs == null)
        {
            CUP_LOGGER.error("sample({}) sig contributions not found", sample.Id);
            return;
        }

        for(Map.Entry<String,Double> entry : sampleSigContribs.entrySet())
        {
            final String sigName = entry.getKey();
            double sampleSigContrib = entry.getValue();

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

            if(!REPORTABLE_SIGS.stream().anyMatch(x -> sigName.equalsIgnoreCase(x)))
                continue;

            for(Map.Entry<String,Map<String,double[]>> cancerContribs : mRefCancerSigContribPercentiles.entrySet())
            {
                final String cancerType = cancerContribs.getKey();
                final double[] refSigPercentiles = cancerContribs.getValue().get(sigName);

                if(refSigPercentiles == null)
                {
                    CUP_LOGGER.error("missing sig({}) data for cancerType({})", sigName, cancerType);
                    return;
                }

                Map<String,Double> cancerResults = Maps.newHashMap();
                double percentile = getPercentile(refSigPercentiles, sampleSigContrib, true);
                cancerResults.put(cancerType, percentile);

                results.add(new SampleResult(sample.Id, SNV_SIG, sigName.toUpperCase(), round(sampleSigContrib), cancerResults));
            }
        }
    }

}
