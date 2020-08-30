package com.hartwig.hmftools.cup.rna;

import static java.lang.Math.max;
import static java.lang.Math.pow;

import static com.hartwig.hmftools.common.sigs.CosineSimilarity.calcCosineSim;
import static com.hartwig.hmftools.common.sigs.SigUtils.loadMatrixDataFile;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CategoryType.CLASSIFIER;
import static com.hartwig.hmftools.cup.common.ClassifierType.RNA_GENE_EXPRESSION;
import static com.hartwig.hmftools.cup.common.CupConstants.RNA_GENE_EXP_CSS_THRESHOLD;
import static com.hartwig.hmftools.cup.common.CupConstants.RNA_GENE_EXP_DIFF_EXPONENT;
import static com.hartwig.hmftools.cup.common.ResultType.LIKELIHOOD;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sigs.SigMatrix;
import com.hartwig.hmftools.cup.SampleAnalyserConfig;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;

public class RnaExpression
{
    private final SampleAnalyserConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private SigMatrix mRefCancerTypeGeneExpression;
    private final List<String> mRefCancerTypes;

    private SigMatrix mSampleRnaExpression;
    private final Map<String,Integer> mSampleIndexMap;

    private boolean mIsValid;

    public RnaExpression(final SampleAnalyserConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mSampleRnaExpression = null;
        mSampleIndexMap = Maps.newHashMap();

        mRefCancerTypeGeneExpression = null;
        mRefCancerTypes = Lists.newArrayList();

        mIsValid = true;

        if(!mConfig.SampleRnaExpFile.isEmpty())
        {
            mSampleRnaExpression = loadMatrixDataFile(mConfig.SampleRnaExpFile, mSampleIndexMap, Lists.newArrayList("GeneId","GeneName"));
            mSampleRnaExpression.cacheTranspose();
        }

        if(!mConfig.RefRnaExpFile.isEmpty())
        {
            mRefCancerTypeGeneExpression = loadMatrixDataFile(mConfig.RefRnaExpFile, mRefCancerTypes, Lists.newArrayList("GeneId","GeneName"));
            mRefCancerTypeGeneExpression.cacheTranspose();
        }
    }

    public boolean isValid() { return mIsValid; }

    public List<SampleResult> processSample(final SampleData sample)
    {
        final List<SampleResult> results = Lists.newArrayList();

        if(mSampleIndexMap.isEmpty())
            return results;

        if(mConfig.runCategory(CLASSIFIER))
        {
            addGeneExpressionCssResults(sample, results);
        }

        return results;
    }

    private void addGeneExpressionCssResults(final SampleData sample, final List<SampleResult> results)
    {
        Integer sampleCountsIndex = mSampleIndexMap.get(sample.Id);

        if(sampleCountsIndex == null)
        {
            // CUP_LOGGER.debug("sample({}) has no RNA data", sample.Id);
            return;
        }

        final double[] sampleGeneTPMs = mSampleRnaExpression.getCol(sampleCountsIndex);

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
                    adjustRefPosFreqCounts(mRefCancerTypeGeneExpression.getCol(i), sampleGeneTPMs) : mRefCancerTypeGeneExpression.getCol(i);

            double css = calcCosineSim(sampleGeneTPMs, refPosFreqs);

            if(css < RNA_GENE_EXP_CSS_THRESHOLD)
                continue;

            double cssWeight = pow(RNA_GENE_EXP_DIFF_EXPONENT, -100 * (1 - css));

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

        results.add(new SampleResult(
                sample.Id, CLASSIFIER, LIKELIHOOD, RNA_GENE_EXPRESSION.toString(), String.format("%.4g", totalCss), cancerCssTotals));
    }

    private double[] adjustRefPosFreqCounts(final double[] refGeneTpmTotals, final double[] sampleGeneTPMs)
    {
        double[] adjustedTPMs = new double[refGeneTpmTotals.length];

        for(int b = 0; b < refGeneTpmTotals.length; ++b)
        {
            adjustedTPMs[b] = max(refGeneTpmTotals[b] - sampleGeneTPMs[b], 0);
        }

        return adjustedTPMs;
    }

}
