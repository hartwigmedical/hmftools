package com.hartwig.hmftools.cup.rna;

import static java.lang.Math.log;
import static java.lang.Math.max;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

import static com.hartwig.hmftools.common.stats.CosineSimilarity.calcCosineSim;
import static com.hartwig.hmftools.common.utils.MatrixUtils.loadMatrixDataFile;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.common.CategoryType.CLASSIFIER;
import static com.hartwig.hmftools.cup.common.CategoryType.GENE_EXP;
import static com.hartwig.hmftools.cup.common.ClassifierType.EXPRESSION_COHORT;
import static com.hartwig.hmftools.cup.common.ClassifierType.EXPRESSION_PAIRWISE;
import static com.hartwig.hmftools.cup.common.CupCalcs.adjustRefCounts;
import static com.hartwig.hmftools.cup.common.CupConstants.RNA_GENE_EXP_CSS_THRESHOLD;
import static com.hartwig.hmftools.cup.common.CupConstants.RNA_GENE_EXP_DIFF_EXPONENT;
import static com.hartwig.hmftools.cup.common.CupConstants.CSS_SIMILARITY_CUTOFF;
import static com.hartwig.hmftools.cup.common.CupConstants.CSS_SIMILARITY_MAX_MATCHES;
import static com.hartwig.hmftools.cup.common.ResultType.LIKELIHOOD;
import static com.hartwig.hmftools.cup.common.SampleData.isKnownCancerType;
import static com.hartwig.hmftools.cup.common.SampleResult.checkIsValidCancerType;
import static com.hartwig.hmftools.cup.common.SampleSimilarity.recordCssSimilarity;
import static com.hartwig.hmftools.cup.rna.RefGeneExpression.loadGeneIdIndices;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.rna.GeneExpressionFile;
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

public class GeneExpressionClassifier implements CuppaClassifier
{
    private final CuppaConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private Matrix mRefCancerTypeGeneExpression;
    private final List<String> mRefCancerTypes;
    private final Map<String,Integer> mRefCancerSampleCounts; // needs to be build from the RNAs if less than the DNA samples available

    private final Map<String,Integer> mGeneIdIndexMap;

    private Matrix mRefSampleGeneExpression;
    private final Map<String,Integer> mRefSampleGeneExpIndexMap;

    private Matrix mSampleRnaExpression;
    private final Map<String,Integer> mSampleIndexMap;

    private final boolean mRunPairwiseCss;
    private final boolean mRunCancerCss;

    private static final String RNA_METHODS = "rna_methods";
    private static final String RNA_METHOD_PAIRWISE_CSS = "pairwise_css";
    private static final String RNA_METHOD_CANCER_CSS = "cancer_css";

    private boolean mIsValid;

    public GeneExpressionClassifier(final CuppaConfig config, final SampleDataCache sampleDataCache, final CommandLine cmd)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mRefSampleGeneExpression = null;
        mRefSampleGeneExpIndexMap = Maps.newHashMap();

        mRefCancerSampleCounts = Maps.newHashMap();

        mRefCancerTypeGeneExpression = null;
        mRefCancerTypes = Lists.newArrayList();
        mGeneIdIndexMap = Maps.newHashMap();

        mSampleRnaExpression = null;
        mSampleIndexMap = Maps.newHashMap();

        final String rnaMethods = cmd.getOptionValue(RNA_METHODS);

        mRunPairwiseCss = rnaMethods == null || rnaMethods.contains(RNA_METHOD_PAIRWISE_CSS);
        mRunCancerCss = rnaMethods != null && rnaMethods.contains(RNA_METHOD_CANCER_CSS);

        mIsValid = true;

        if(mSampleDataCache.isMultiSample() && mConfig.SampleGeneExpFile.isEmpty())
            return;

        if(mRunPairwiseCss && mConfig.RefGeneExpSampleFile.isEmpty())
            return;

        if(mRunCancerCss && mConfig.RefGeneExpCancerFile.isEmpty())
            return;

        final List<String> ignoreFields = Lists.newArrayList("GeneId", "GeneName");

        if(mRunPairwiseCss)
        {
            loadGeneIdIndices(mConfig.RefGeneExpSampleFile, mGeneIdIndexMap);
            mRefSampleGeneExpression = loadMatrixDataFile(mConfig.RefGeneExpSampleFile, mRefSampleGeneExpIndexMap, ignoreFields);

            if(mRefSampleGeneExpression ==  null)
            {
                mIsValid = false;
                return;
            }

            mRefSampleGeneExpression.cacheTranspose();
        }

        if(mRunCancerCss)
        {
            if(mGeneIdIndexMap.isEmpty())
                loadGeneIdIndices(mConfig.RefGeneExpCancerFile, mGeneIdIndexMap);

            mRefCancerTypeGeneExpression = loadMatrixDataFile(mConfig.RefGeneExpCancerFile, mRefCancerTypes, ignoreFields);

            if(mRefCancerTypeGeneExpression ==  null)
            {
                mIsValid = false;
                return;
            }

            mRefCancerTypeGeneExpression.cacheTranspose();
        }

        buildCancerSampleCounts();

        if(mConfig.SampleGeneExpFile.equals(mConfig.RefGeneExpSampleFile))
        {
            mSampleRnaExpression = mRefSampleGeneExpression;
            mSampleIndexMap.putAll(mRefSampleGeneExpIndexMap);
        }
        else
        {
            if(mSampleDataCache.isSingleSample())
            {
                loadSampleGeneExpressionData(mSampleDataCache.SampleIds.get(0));
            }
            else
            {
                mSampleRnaExpression = loadMatrixDataFile(mConfig.SampleGeneExpFile, mSampleIndexMap, ignoreFields);

                if(mSampleRnaExpression == null)
                {
                    mIsValid = false;
                    return;
                }

                mSampleRnaExpression.cacheTranspose();
            }
        }
    }

    private void buildCancerSampleCounts()
    {
        if(mRefSampleGeneExpIndexMap.isEmpty())
        {
            for(Map.Entry<String,List<SampleData>> entry : mSampleDataCache.RefCancerSampleData.entrySet())
            {
                mRefCancerSampleCounts.put(entry.getKey(), entry.getValue().size());
            }
        }
        else
        {
            for(final String refSampleId : mRefSampleGeneExpIndexMap.keySet())
            {
                final String refCancerType = mSampleDataCache.RefSampleCancerTypeMap.get(refSampleId);

                if(refCancerType == null)
                    continue;

                Integer sampleCount = mRefCancerSampleCounts.get(refCancerType);
                if(sampleCount == null)
                    mRefCancerSampleCounts.put(refCancerType, 1);
                else
                    mRefCancerSampleCounts.put(refCancerType, sampleCount + 1);
            }
        }

        for(Map.Entry<String,Integer> entry : mRefCancerSampleCounts.entrySet())
        {
            CUP_LOGGER.debug("RNA ref cancer({}) samples({})", entry.getKey(), entry.getValue());
        }
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(RNA_METHODS, true, "Types of RNA gene expression methods");
    }

    public CategoryType categoryType() { return GENE_EXP; }
    public boolean isValid() { return mIsValid; }
    public void close() {}

    public void processSample(final SampleData sample, final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        if(mSampleIndexMap.isEmpty())
            return;

        if(mSampleDataCache.isMultiSample() && !sample.hasRna())
            return;

        Integer sampleCountsIndex = mSampleIndexMap.get(sample.Id);

        if(sampleCountsIndex == null)
        {
            CUP_LOGGER.warn("sample({}) gene expression matrix data not found", sample.Id);
            return;
        }

        final double[] sampleGeneTPMs = mSampleRnaExpression.getCol(sampleCountsIndex);

        if(mRunCancerCss)
            addCancerCssResults(sample, sampleGeneTPMs, results);

        if(mRunPairwiseCss)
            addSampleCssResults(sample, sampleGeneTPMs, results, similarities);
    }

    private void addCancerCssResults(final SampleData sample, final double[] sampleGeneTPMs, final List<SampleResult> results)
    {
        int refCancerCount = mRefCancerTypeGeneExpression.Cols;

        final Map<String,Double> cancerCssTotals = Maps.newHashMap();

        for(int i = 0; i < refCancerCount; ++i)
        {
            final String refCancerType = mRefCancerTypes.get(i);

            if(!isKnownCancerType(refCancerType))
                continue;

            if(!checkIsValidCancerType(sample, refCancerType, cancerCssTotals))
                continue;

            boolean matchesCancerType = sample.cancerType().equals(refCancerType);

            final double[] refPosFreqs = sample.isRefSample() && matchesCancerType ?
                    adjustRefCounts(mRefCancerTypeGeneExpression.getCol(i), sampleGeneTPMs, 1) : mRefCancerTypeGeneExpression.getCol(i);

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
                sample.Id, CLASSIFIER, LIKELIHOOD, EXPRESSION_COHORT.toString(), String.format("%.4g", totalCss), cancerCssTotals));
    }

    private void addSampleCssResults(
            final SampleData sample, final double[] sampleTPMs, final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        final Map<String,Double> cancerCssTotals = Maps.newHashMap();

        final List<SampleSimilarity> topMatches = Lists.newArrayList();

        for(Map.Entry<String,Integer> entry : mRefSampleGeneExpIndexMap.entrySet())
        {
            final String refSampleId = entry.getKey();

            if(refSampleId.equals(sample.Id))
                continue;

            final String refCancerType = mSampleDataCache.RefSampleCancerTypeMap.get(refSampleId);

            if(refCancerType == null)
                continue;

            if(!checkIsValidCancerType(sample, refCancerType, cancerCssTotals))
                continue;

            int refSampleCountsIndex = entry.getValue();
            final double[] otherSampleTPMs = mRefSampleGeneExpression.getCol(refSampleCountsIndex);

            double css = calcCosineSim(sampleTPMs, otherSampleTPMs);

            if(css < RNA_GENE_EXP_CSS_THRESHOLD)
                continue;

            if(mConfig.WriteSimilarities)
            {
                recordCssSimilarity(
                        topMatches, sample.Id, refSampleId, css, EXPRESSION_PAIRWISE.toString(),
                        CSS_SIMILARITY_MAX_MATCHES, CSS_SIMILARITY_CUTOFF);
            }

            if(!isKnownCancerType(refCancerType))
                continue;

            double cssWeight = pow(RNA_GENE_EXP_DIFF_EXPONENT, -100 * (1 - css));

            int cancerSampleCount = mRefCancerSampleCounts.get(refCancerType);
            double weightedCss = css * cssWeight / sqrt(cancerSampleCount);

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
                sample.Id, CLASSIFIER, LIKELIHOOD, EXPRESSION_PAIRWISE.toString(), String.format("%.4g", totalCss), cancerCssTotals));

        similarities.addAll(topMatches);
    }

    private void loadSampleGeneExpressionData(final String sampleId)
    {
        final String filename = GeneExpressionFile.generateFilename(mConfig.SampleDataDir, sampleId);

        if(!Files.exists(Paths.get(filename)))
            return;

        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());
            String header = fileData.get(0);
            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, ",");
            fileData.remove(0);

            mSampleRnaExpression = new Matrix(mGeneIdIndexMap.size(), 1);
            mSampleIndexMap.put(mSampleDataCache.SampleIds.get(0), 0);

            // GeneId,GeneName, etc AdjTPM

            int geneIdCol = fieldsIndexMap.get("GeneId");
            int adjTPM = fieldsIndexMap.get("AdjTPM");

            for(String line : fileData)
            {
                final String[] items = line.split(DATA_DELIM, -1);

                String geneId = items[geneIdCol];
                double adjTpm = Double.parseDouble(items[adjTPM]);
                Integer geneIdIndex = mGeneIdIndexMap.get(geneId);

                if(geneIdIndex == null)
                {
                    CUP_LOGGER.error("unknown geneId({}) in sample file({})", geneId, filename);
                    return;
                }

                double logTpm = log(adjTpm + 1);

                mSampleRnaExpression.set(geneIdIndex, 0, logTpm);
            }

            mSampleRnaExpression.cacheTranspose();
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read RNA sample gene data file({}): {}", filename, e.toString());
            return;
        }
    }

}
