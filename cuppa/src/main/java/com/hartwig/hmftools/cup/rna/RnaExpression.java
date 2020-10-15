package com.hartwig.hmftools.cup.rna;

import static java.lang.Math.exp;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

import static com.hartwig.hmftools.common.sigs.CosineSimilarity.calcCosineSim;
import static com.hartwig.hmftools.common.sigs.SigUtils.loadMatrixDataFile;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CategoryType.CLASSIFIER;
import static com.hartwig.hmftools.cup.common.CategoryType.GENE_EXP;
import static com.hartwig.hmftools.cup.common.CategoryType.SV;
import static com.hartwig.hmftools.cup.common.ClassifierType.RNA_GENE_EXPRESSION;
import static com.hartwig.hmftools.cup.common.ClassifierType.RNA_GENE_EXPRESSION_PAIRWISE;
import static com.hartwig.hmftools.cup.common.CupCalcs.calcPercentilePrevalence;
import static com.hartwig.hmftools.cup.common.CupCalcs.convertToPercentages;
import static com.hartwig.hmftools.cup.common.CupConstants.RNA_GENE_EXP_CSS_THRESHOLD;
import static com.hartwig.hmftools.cup.common.CupConstants.RNA_GENE_EXP_DIFF_EXPONENT;
import static com.hartwig.hmftools.cup.common.ResultType.LIKELIHOOD;
import static com.hartwig.hmftools.cup.common.ResultType.PERCENTILE;
import static com.hartwig.hmftools.cup.rna.RefRnaExpression.loadRefPercentileData;
import static com.hartwig.hmftools.cup.rna.RefRnaExpression.populateGeneIdList;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sigs.SigMatrix;
import com.hartwig.hmftools.cup.SampleAnalyserConfig;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;
import com.hartwig.hmftools.cup.svs.SvData;
import com.hartwig.hmftools.cup.svs.SvDataType;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class RnaExpression
{
    private final SampleAnalyserConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private SigMatrix mRefCancerTypeGeneExpression;
    private final List<String> mRefCancerTypes;

    private final Map<String,Map<String,double[]>> mRefGeneCancerPercentiles;
    private final Map<String,String> mGeneIdNameMap;
    private final List<String> mGeneIdList;

    private SigMatrix mSampleRnaExpression;
    private final Map<String,Integer> mSampleIndexMap;

    private final boolean mIncludePairwiseCss;
    private static final String INCLUDE_PAIRWISE_CSS = "rna_pairwise_css";

    private final boolean mExcludeCancerCss;
    private static final String EXCLUDE_CANCER_CSS = "rna_exclude_cancer_css";

    private boolean mIsValid;

    public RnaExpression(final SampleAnalyserConfig config, final SampleDataCache sampleDataCache, final CommandLine cmd)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mSampleRnaExpression = null;
        mSampleIndexMap = Maps.newHashMap();

        mRefCancerTypeGeneExpression = null;
        mRefCancerTypes = Lists.newArrayList();
        mGeneIdList = Lists.newArrayList();

        mRefGeneCancerPercentiles = Maps.newHashMap();
        mGeneIdNameMap = Maps.newHashMap();

        mIncludePairwiseCss = cmd.hasOption(INCLUDE_PAIRWISE_CSS);
        mExcludeCancerCss = cmd.hasOption(EXCLUDE_CANCER_CSS);

        mIsValid = true;

        if(mConfig.SampleRnaExpFile.isEmpty() && mConfig.RefRnaCancerExpFile.isEmpty())
            return;

        mSampleRnaExpression = loadMatrixDataFile(mConfig.SampleRnaExpFile, mSampleIndexMap, Lists.newArrayList("GeneId","GeneName"));
        mRefCancerTypeGeneExpression = loadMatrixDataFile(mConfig.RefRnaCancerExpFile, mRefCancerTypes, Lists.newArrayList("GeneId","GeneName"));

        if(!mConfig.RefRnaGeneCancerPercFile.isEmpty())
        {
            mIsValid &= loadRefPercentileData(mConfig.RefRnaGeneCancerPercFile, mRefGeneCancerPercentiles, mGeneIdNameMap);
            populateGeneIdList(mConfig.RefRnaCancerExpFile, mGeneIdList);
        }

        if(mSampleRnaExpression == null || mRefCancerTypeGeneExpression ==  null)
        {
            mIsValid = false;
            return;
        }

        mSampleRnaExpression.cacheTranspose();
        mRefCancerTypeGeneExpression.cacheTranspose();
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(INCLUDE_PAIRWISE_CSS, false, "Run RNA gene expression sample pairwise CSS");
        options.addOption(EXCLUDE_CANCER_CSS, false, "Exclude RNA gene expression cancer CSS");
    }

    public boolean isValid() { return mIsValid; }

    public List<SampleResult> processSample(final SampleData sample)
    {
        final List<SampleResult> results = Lists.newArrayList();

        if(mSampleIndexMap.isEmpty())
            return results;

        Integer sampleCountsIndex = mSampleIndexMap.get(sample.Id);

        if(sampleCountsIndex == null)
            return results;

        final double[] sampleGeneTPMs = mSampleRnaExpression.getCol(sampleCountsIndex);

        if(!mExcludeCancerCss)
            addCancerCssResults(sample, sampleGeneTPMs, results);

        if(mIncludePairwiseCss && mSampleRnaExpression.Cols == mSampleDataCache.RefSampleCancerTypeMap.size())
            addSampleCssResults(sample, sampleGeneTPMs, results);

        if(!mRefGeneCancerPercentiles.isEmpty())
            addPrevalenceResults(sample, sampleGeneTPMs, results);

        return results;
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
                sample.Id, CLASSIFIER, LIKELIHOOD, RNA_GENE_EXPRESSION.toString(), String.format("%.4g", totalCss), cancerCssTotals));
    }

    private void addSampleCssResults(final SampleData sample, final double[] sampleTPMs, final List<SampleResult> results)
    {
        final Map<String,Double> cancerCssTotals = Maps.newHashMap();

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

            double cssWeight = pow(RNA_GENE_EXP_DIFF_EXPONENT, -100 * (1 - css));

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
            double prob = totalCss > 0 ? entry.getValue() / totalCss : 0;
            cancerCssTotals.put(entry.getKey(), prob);
        }

        results.add(new SampleResult(
                sample.Id, CLASSIFIER, LIKELIHOOD, RNA_GENE_EXPRESSION_PAIRWISE.toString(), String.format("%.4g", totalCss), cancerCssTotals));
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
        int cancerSampleCount = sample.isRefSample() ? mSampleDataCache.RefCancerSampleData.get(sample.CancerType).size() : 0;

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
