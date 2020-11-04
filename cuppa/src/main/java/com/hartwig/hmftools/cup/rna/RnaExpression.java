package com.hartwig.hmftools.cup.rna;

import static java.lang.Math.exp;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

import static com.hartwig.hmftools.common.sigs.CosineSimilarity.calcCosineSim;
import static com.hartwig.hmftools.common.sigs.SigUtils.loadMatrixDataFile;
import static com.hartwig.hmftools.cup.common.CategoryType.CLASSIFIER;
import static com.hartwig.hmftools.cup.common.CategoryType.FEATURE;
import static com.hartwig.hmftools.cup.common.CategoryType.GENE_EXP;
import static com.hartwig.hmftools.cup.common.ClassifierType.GENE_EXPRESSION_CANCER;
import static com.hartwig.hmftools.cup.common.ClassifierType.GENE_EXPRESSION_PAIRWISE;
import static com.hartwig.hmftools.cup.common.CupCalcs.calcPercentilePrevalence;
import static com.hartwig.hmftools.cup.common.CupCalcs.convertToPercentages;
import static com.hartwig.hmftools.cup.common.CupConstants.RNA_GENE_EXP_CSS_THRESHOLD;
import static com.hartwig.hmftools.cup.common.CupConstants.RNA_GENE_EXP_DIFF_EXPONENT;
import static com.hartwig.hmftools.cup.common.CupConstants.CSS_SIMILARITY_CUTOFF;
import static com.hartwig.hmftools.cup.common.CupConstants.CSS_SIMILARITY_MAX_MATCHES;
import static com.hartwig.hmftools.cup.common.ResultType.LIKELIHOOD;
import static com.hartwig.hmftools.cup.common.SampleSimilarity.recordCssSimilarity;
import static com.hartwig.hmftools.cup.rna.RefRnaExpression.loadRefPercentileData;
import static com.hartwig.hmftools.cup.rna.RefRnaExpression.populateGeneIdList;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sigs.SigMatrix;
import com.hartwig.hmftools.cup.CuppaConfig;
import com.hartwig.hmftools.cup.common.CategoryType;
import com.hartwig.hmftools.cup.common.CuppaClassifier;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;
import com.hartwig.hmftools.cup.common.SampleSimilarity;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class RnaExpression implements CuppaClassifier
{
    private final CuppaConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private SigMatrix mRefCancerTypeGeneExpression;
    private final List<String> mRefCancerTypes;

    private final Map<String,Map<String,double[]>> mRefGeneCancerPercentiles;
    private final Map<String,String> mGeneIdNameMap;
    private final List<String> mGeneIdList;

    private SigMatrix mRefSampleGeneExpression;
    private final Map<String,Integer> mRefSampleGeneExpIndexMap;

    private SigMatrix mSampleRnaExpression;
    private final Map<String,Integer> mSampleIndexMap;

    private final boolean mRunPairwiseCss;
    private final boolean mRunCancerCss;
    private final boolean mRunGenePrevalence;

    private static final String RNA_METHODS = "rna_methods";
    private static final String RNA_METHOD_PAIRWISE_CSS = "pairwise_css";
    private static final String RNA_METHOD_CANCER_CSS = "cancer_css";
    private static final String RNA_METHOD_GENE_PREV = "gene_prevalence";

    private boolean mIsValid;

    public RnaExpression(final CuppaConfig config, final SampleDataCache sampleDataCache, final CommandLine cmd)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mRefSampleGeneExpression = null;
        mRefSampleGeneExpIndexMap = Maps.newHashMap();

        mRefCancerTypeGeneExpression = null;
        mRefCancerTypes = Lists.newArrayList();
        mGeneIdList = Lists.newArrayList();

        mRefGeneCancerPercentiles = Maps.newHashMap();
        mGeneIdNameMap = Maps.newHashMap();

        mSampleRnaExpression = null;
        mSampleIndexMap = Maps.newHashMap();

        final String rnaMethods = cmd.getOptionValue(RNA_METHODS);

        mRunPairwiseCss = rnaMethods != null && rnaMethods.contains(RNA_METHOD_PAIRWISE_CSS);
        mRunCancerCss = rnaMethods == null || rnaMethods.contains(RNA_METHOD_CANCER_CSS);
        mRunGenePrevalence = rnaMethods != null && rnaMethods.contains(RNA_METHOD_GENE_PREV);

        mIsValid = true;

        if(mConfig.SampleRnaExpFile.isEmpty())
            return;

        if(mRunPairwiseCss && mConfig.RefGeneExpSampleFile.isEmpty())
            return;

        if(mRunCancerCss && mConfig.RefGeneExpCancerFile.isEmpty())
            return;

        if(mRunGenePrevalence && mConfig.RefGeneExpPercFile.isEmpty())
            return;

        if(!mConfig.RefGeneExpSampleFile.isEmpty())
        {
            mRefSampleGeneExpression =
                    loadMatrixDataFile(mConfig.RefGeneExpSampleFile, mRefSampleGeneExpIndexMap, Lists.newArrayList("GeneId", "GeneName"));

            if(mRefSampleGeneExpression ==  null)
            {
                mIsValid = false;
                return;
            }

            mRefSampleGeneExpression.cacheTranspose();
        }

        if(!mConfig.RefGeneExpCancerFile.isEmpty())
        {
            mRefCancerTypeGeneExpression =
                    loadMatrixDataFile(mConfig.RefGeneExpCancerFile, mRefCancerTypes, Lists.newArrayList("GeneId", "GeneName"));

            if(mRefCancerTypeGeneExpression ==  null)
            {
                mIsValid = false;
                return;
            }

            mRefCancerTypeGeneExpression.cacheTranspose();
        }

        if(!mConfig.RefGeneExpPercFile.isEmpty())
        {
            mIsValid &= loadRefPercentileData(mConfig.RefGeneExpPercFile, mRefGeneCancerPercentiles, mGeneIdNameMap);
            populateGeneIdList(mConfig.RefGeneExpCancerFile, mGeneIdList);
        }

        if(mConfig.SampleRnaExpFile.equals(mConfig.RefGeneExpSampleFile))
        {
            mSampleRnaExpression = mRefSampleGeneExpression;
            mSampleIndexMap.putAll(mRefSampleGeneExpIndexMap);
        }
        else
        {
            mSampleRnaExpression = loadMatrixDataFile(mConfig.SampleRnaExpFile, mSampleIndexMap, Lists.newArrayList("GeneId", "GeneName"));

            if(mSampleRnaExpression == null || mRefCancerTypeGeneExpression ==  null)
            {
                mIsValid = false;
                return;
            }

            mSampleRnaExpression.cacheTranspose();
        }
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(RNA_METHODS, true, "Types of RNA gene expression methods");
    }

    public CategoryType categoryType() { return GENE_EXP; }
    public boolean isValid() { return mIsValid; }

    public void processSample(final SampleData sample, final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        if(mSampleIndexMap.isEmpty())
            return;

        Integer sampleCountsIndex = mSampleIndexMap.get(sample.Id);

        if(sampleCountsIndex == null)
            return;

        final double[] sampleGeneTPMs = mSampleRnaExpression.getCol(sampleCountsIndex);

        if(mRunCancerCss)
            addCancerCssResults(sample, sampleGeneTPMs, results);

        if(mRunPairwiseCss && mSampleRnaExpression.Cols == mSampleDataCache.RefSampleCancerTypeMap.size())
            addSampleCssResults(sample, sampleGeneTPMs, results, similarities);

        if(mRunGenePrevalence && !mRefGeneCancerPercentiles.isEmpty())
            addPrevalenceResults(sample, sampleGeneTPMs, results);
    }

    private void addCancerCssResults(final SampleData sample, final double[] sampleGeneTPMs, final List<SampleResult> results)
    {
        int refCancerCount = mRefCancerTypeGeneExpression.Cols;

        final Map<String,Double> cancerCssTotals = Maps.newHashMap();

        for(int i = 0; i < refCancerCount; ++i)
        {
            final String refCancerType = mRefCancerTypes.get(i);

            if(!sample.isCandidateCancerType(refCancerType))
            {
                cancerCssTotals.put(refCancerType, 0.0);
                continue;
            }

            boolean matchesCancerType = sample.CancerType.equals(refCancerType);

            final double[] refPosFreqs = sample.isRefSample() && matchesCancerType ?
                    adjustRefTpmTotals(mRefCancerTypeGeneExpression.getCol(i), sampleGeneTPMs) : mRefCancerTypeGeneExpression.getCol(i);

            double css = calcCosineSim(sampleGeneTPMs, refPosFreqs);

            if(css < RNA_GENE_EXP_CSS_THRESHOLD)
                continue;

            double cssWeight = pow(RNA_GENE_EXP_DIFF_EXPONENT, -100 * (1 - css));

            double weightedCss = css * cssWeight;
            cancerCssTotals.put(refCancerType, weightedCss);
        }

        double totalCss = cancerCssTotals.values().stream().mapToDouble(x -> x).sum();

        for(Map.Entry<String,Double> entry : cancerCssTotals.entrySet())
        {
            cancerCssTotals.put(entry.getKey(), entry.getValue() / totalCss);
        }

        results.add(new SampleResult(
                sample.Id, CLASSIFIER, LIKELIHOOD, GENE_EXPRESSION_CANCER.toString(), String.format("%.4g", totalCss), cancerCssTotals));
    }

    private void addSampleCssResults(
            final SampleData sample, final double[] sampleTPMs, final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        final Map<String,Double> cancerCssTotals = Maps.newHashMap();

        final List<SampleSimilarity> topMatches = Lists.newArrayList();

        for(Map.Entry<String,Integer> entry : mSampleIndexMap.entrySet())
        {
            final String refSampleId = entry.getKey();

            if(refSampleId.equals(sample.Id))
                continue;

            final String refCancerType = mSampleDataCache.RefSampleCancerTypeMap.get(refSampleId);

            if(refCancerType == null)
                continue;

            if(!sample.isCandidateCancerType(refCancerType))
            {
                cancerCssTotals.put(refCancerType, 0.0);
                continue;
            }

            int refSampleCountsIndex = entry.getValue();
            final double[] otherSampleTPMs = mSampleRnaExpression.getCol(refSampleCountsIndex);

            double css = calcCosineSim(sampleTPMs, otherSampleTPMs);

            if(css < RNA_GENE_EXP_CSS_THRESHOLD)
                continue;

            if(mConfig.WriteSimilarities)
            {
                recordCssSimilarity(
                        topMatches, sample.Id, refSampleId, css, GENE_EXPRESSION_PAIRWISE.toString(),
                        CSS_SIMILARITY_MAX_MATCHES, CSS_SIMILARITY_CUTOFF);
            }

            double cssWeight = pow(RNA_GENE_EXP_DIFF_EXPONENT, -100 * (1 - css));

            int cancerTypeCount = mSampleDataCache.getCancerSampleCount(refCancerType);
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
            double prob = totalCss > 0 ? entry.getValue() / totalCss : 0;
            cancerCssTotals.put(entry.getKey(), prob);
        }

        results.add(new SampleResult(
                sample.Id, CLASSIFIER, LIKELIHOOD, GENE_EXPRESSION_PAIRWISE.toString(), String.format("%.4g", totalCss), cancerCssTotals));

        similarities.addAll(topMatches);
    }

    private double[] adjustRefTpmTotals(final double[] refGeneTpmTotals, final double[] sampleGeneTPMs)
    {
        double[] adjustedTPMs = new double[refGeneTpmTotals.length];

        for(int b = 0; b < refGeneTpmTotals.length; ++b)
        {
            adjustedTPMs[b] = max(refGeneTpmTotals[b] - sampleGeneTPMs[b], 0);
        }

        return adjustedTPMs;
    }

    private static final double MAX_PERC_THRESHOLD = 0.50;
    private static final double MIN_TPM_THRESHOLD = 0.1;

    private void addPrevalenceResults(final SampleData sample, final double[] sampleGeneTPMs, final List<SampleResult> results)
    {
        int cancerTypeCount = mSampleDataCache.RefCancerSampleData.size();
        int cancerSampleCount = sample.isRefSample() ? mSampleDataCache.getCancerSampleCount(sample.CancerType) : 0;

        final Map<String,Double> summaryCancerPrevs = Maps.newHashMap();

        for(int i = 0; i < mGeneIdList.size(); ++i)
        {
            final String geneId = mGeneIdList.get(i);
            final String geneName = mGeneIdNameMap.get(geneId);

            if(geneName == null)
                continue;

            double sampleLogTpm = sampleGeneTPMs[i];
            double sampleTpm = exp(sampleLogTpm) - 1; // revert to TPM from log(TPM+1)

            if(sampleTpm < MIN_TPM_THRESHOLD)
                continue;

            final Map<String,double[]> cancerPercentiles = mRefGeneCancerPercentiles.get(geneId);

            final Map<String,Double> cancerPrevs = calcPercentilePrevalence(
                    sample.CancerType, cancerSampleCount,  cancerTypeCount, cancerPercentiles, sampleTpm, false);


            /*
            double maxPercentage = cancerPrevs.values().stream().mapToDouble(x -> x).max().orElse(0);
            final String dataType = String.format("%s:%s", geneId, geneName);

            if(maxPercentage >= MAX_PERC_THRESHOLD)
            {
                results.add(new SampleResult(sample.Id, GENE_EXP, PERCENTILE, dataType, String.format("%.3g", sampleTpm), cancerPrevs));
            }
            */

            for(Map.Entry<String,Double> entry : cancerPrevs.entrySet())
            {
                Double prev = summaryCancerPrevs.get(entry.getKey());
                if(prev == null)
                    summaryCancerPrevs.put(entry.getKey(), entry.getValue());
                else
                    summaryCancerPrevs.put(entry.getKey(), prev * entry.getValue());
            }

            double maxProbability = summaryCancerPrevs.values().stream().mapToDouble(x -> x).max().orElse(0);

            if(maxProbability <= 0)
                break;

            if(maxProbability < 1e-20)
                convertToPercentages(summaryCancerPrevs); // prevent values getting too small and dropping out
        }

        // form a likelihood value from all gene entries
        double totalPrevalence = summaryCancerPrevs.values().stream().mapToDouble(x -> x).sum();

        if(totalPrevalence > 0)
        {
            convertToPercentages(summaryCancerPrevs);
            results.add(new SampleResult(sample.Id, GENE_EXP, LIKELIHOOD, "GENE_EXP_LIKELIHOOD", "", summaryCancerPrevs));
        }
    }
}
