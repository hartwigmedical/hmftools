package com.hartwig.hmftools.cup.rna;

import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.stats.CosineSimilarity.calcCosineSim;
import static com.hartwig.hmftools.common.utils.MatrixFile.loadMatrixDataFile;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CategoryType.GENE_EXP;
import static com.hartwig.hmftools.cup.common.ClassifierType.EXPRESSION_COHORT;
import static com.hartwig.hmftools.cup.common.ClassifierType.EXPRESSION_PAIRWISE;
import static com.hartwig.hmftools.cup.common.CupCalcs.adjustRefCounts;
import static com.hartwig.hmftools.cup.common.CupConstants.GENE_EXP_CSS_THRESHOLD;
import static com.hartwig.hmftools.cup.common.CupConstants.GENE_EXP_DIFF_EXPONENT;
import static com.hartwig.hmftools.cup.common.CupConstants.CSS_SIMILARITY_CUTOFF;
import static com.hartwig.hmftools.cup.common.CupConstants.CSS_SIMILARITY_MAX_MATCHES;
import static com.hartwig.hmftools.cup.common.ResultType.CLASSIFIER;
import static com.hartwig.hmftools.cup.common.SampleData.isKnownCancerType;
import static com.hartwig.hmftools.cup.common.SampleResult.checkIsValidCancerType;
import static com.hartwig.hmftools.cup.common.SampleSimilarity.recordCssSimilarity;
import static com.hartwig.hmftools.cup.rna.GeneExpressionDataLoader.GENE_EXP_IGNORE_FIELDS;
import static com.hartwig.hmftools.cup.rna.GeneExpressionDataLoader.loadGeneExpressionMatrix;
import static com.hartwig.hmftools.cup.rna.GeneExpressionDataLoader.loadGeneIdIndices;
import static com.hartwig.hmftools.cup.rna.GeneExpressionDataLoader.loadSampleGeneExpressionFile;
import static com.hartwig.hmftools.cup.rna.GeneExpressionDataLoader.loadSampleGeneExpressionMatrix;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.rna.GeneExpressionFile;
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

public class GeneExpressionClassifier implements CuppaClassifier
{
    private final CuppaConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private Matrix mRefCancerTypeGeneExpression;
    private final List<String> mRefCancerTypes;
    private final Map<String,Integer> mRefCancerSampleCounts; // needs to be build from the RNAs if less than the DNA samples available

    private final Map<String,Integer> mGeneIdIndexMap; // to ensure genes are ordered consistently in ref & sample matrices

    private Matrix mRefSampleGeneExpression;
    private final Map<String,Integer> mRefSampleGeneExpIndexMap;

    private Matrix mSampleGeneExpression;
    private final Map<String,Integer> mSampleIndexMap;

    private final boolean mRunPairwiseCss;
    private final boolean mRunCancerCss;
    private final double mCssExponent;
    private final boolean mMatchReadLength;

    private static final String RNA_METHODS = "gene_exp_rna_methods";
    private static final String CSS_METHOD_PAIRWISE = "pairwise_css";
    private static final String CSS_METHOD_CANCER = "cancer_css";
    private static final String CSS_EXPONENT = "gene_exp_css_exp";
    private static final String MATCH_READ_LENGTH = "gene_exp_match_read_length";

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

        mSampleGeneExpression = null;
        mSampleIndexMap = Maps.newHashMap();

        final String rnaMethods = cmd.getOptionValue(RNA_METHODS);

        mRunPairwiseCss = rnaMethods == null || rnaMethods.contains(CSS_METHOD_PAIRWISE);
        mRunCancerCss = rnaMethods != null && rnaMethods.contains(CSS_METHOD_CANCER);
        mCssExponent = Double.parseDouble(cmd.getOptionValue(CSS_EXPONENT, String.valueOf(GENE_EXP_DIFF_EXPONENT)));
        mMatchReadLength = cmd.hasOption(MATCH_READ_LENGTH);

        mIsValid = true;

        if(mRunPairwiseCss && mConfig.RefGeneExpSampleFile.isEmpty())
            return;

        if(mRunCancerCss && mConfig.RefGeneExpCancerFile.isEmpty())
            return;

        final List<String> sampleNames = Lists.newArrayList();
        final List<String> geneNames = Lists.newArrayList();
        final List<String> geneIds = Lists.newArrayList();

        mRefSampleGeneExpression = loadGeneExpressionMatrix(
                mConfig.RefGeneExpSampleFile, mRefSampleGeneExpIndexMap, sampleNames, geneIds, geneNames);

        for(int i = 0; i < geneIds.size(); ++i)
        {
            mGeneIdIndexMap.put(geneIds.get(i), i);
        }

        if(mRefSampleGeneExpression ==  null)
        {
            mIsValid = false;
            return;
        }

        if(mRunCancerCss)
        {
            // TODO: could make a stand-alone method to load in column headers and geneIds, and use this for both cancer and sample data
            if(mGeneIdIndexMap.isEmpty())
                loadGeneIdIndices(mConfig.RefGeneExpCancerFile, mGeneIdIndexMap);

            mRefCancerTypeGeneExpression = loadMatrixDataFile(mConfig.RefGeneExpCancerFile, mRefCancerTypes, GENE_EXP_IGNORE_FIELDS, true);

            if(mRefCancerTypeGeneExpression ==  null)
            {
                mIsValid = false;
                return;
            }
        }

        buildCancerSampleCounts();

        if(mConfig.TestRefData)
        {
            CUP_LOGGER.debug("re-using ref sample gene-expression matrix data for samples");
            mSampleGeneExpression = mRefSampleGeneExpression;
            mSampleIndexMap.putAll(mRefSampleGeneExpIndexMap);
        }
        else
        {
            for(SampleData sample : mSampleDataCache.SampleDataList)
            {
                final String isofoxDir = mConfig.getIsofoxDataDir(sample.Id);
                final String filename = GeneExpressionFile.generateFilename(isofoxDir, sample.Id);
                CUP_LOGGER.debug("loading sample gene-expression data file({})", filename);

                mSampleIndexMap.put(sample.Id, 0);
                mSampleGeneExpression = loadSampleGeneExpressionFile(filename, mGeneIdIndexMap);
            }

            if(mSampleGeneExpression == null)
            {
                mIsValid = false;
                return;
            }
        }

        // apply any specified noise
        if(mConfig.NoiseAdjustments.makeNoiseAdjustment(EXPRESSION_PAIRWISE))
        {
            final double[] noiseAdjustments = mConfig.NoiseAdjustments.getNoiseData(EXPRESSION_PAIRWISE);
            int noiseAllocation = mConfig.NoiseAdjustments.getNoiseAllocation(EXPRESSION_PAIRWISE);

            CUP_LOGGER.debug("appying noise({}) to gene expression", noiseAllocation);

            NoiseRefCache.applyNoise(mRefSampleGeneExpression, noiseAdjustments, noiseAllocation);

            if(mSampleGeneExpression != mRefSampleGeneExpression)
                NoiseRefCache.applyNoise(mSampleGeneExpression, noiseAdjustments, noiseAllocation);
        }
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(RNA_METHODS, true, "Types of RNA gene expression methods");
        options.addOption(CSS_EXPONENT, true, "Gene expression CSS exponent");
        options.addOption(MATCH_READ_LENGTH, false, "Gene expression pairwise only amongst matching read-length samples");
    }

    public CategoryType categoryType() { return GENE_EXP; }
    public boolean isValid() { return mIsValid; }
    public void close() {}

    private String formCancerType(final SampleData sample)
    {
        if(mMatchReadLength)
            return format("%s_%d", sample.cancerType(), sample.rnaReadLength());
        else
            return sample.cancerType();
    }

    private void buildCancerSampleCounts()
    {
        for(Map.Entry<String,List<SampleData>> refCancerEntry : mSampleDataCache.RefCancerSampleData.entrySet())
        {
            for(SampleData refSample : refCancerEntry.getValue())
            {
                if(!mRefSampleGeneExpIndexMap.containsKey(refSample.Id))
                    continue;

                final String cancerType = formCancerType(refSample);

                Integer sampleCount = mRefCancerSampleCounts.get(cancerType);
                if(sampleCount == null)
                    mRefCancerSampleCounts.put(cancerType, 1);
                else
                    mRefCancerSampleCounts.put(cancerType, sampleCount + 1);
            }
        }

        for(Map.Entry<String,Integer> entry : mRefCancerSampleCounts.entrySet())
        {
            CUP_LOGGER.debug("RNA ref cancer({}) samples({})", entry.getKey(), entry.getValue());
        }
    }

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

        final double[] sampleGeneTPMs = mSampleGeneExpression.getRow(sampleCountsIndex);

        if(mRunCancerCss)
            addCancerCssResults(sample, sampleGeneTPMs, results);

        if(mRunPairwiseCss)
            addSampleCssResults(sample, sampleGeneTPMs, results, similarities);
    }

    private void addSampleCssResults(
            final SampleData sample, final double[] sampleTPMs, final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        final Map<String,Double> cancerCssTotals = Maps.newHashMap();

        final List<SampleSimilarity> topMatches = Lists.newArrayList();

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

                if(mMatchReadLength && refSample.rnaReadLength() != sample.rnaReadLength())
                    continue;

                Integer refSampleIndex = mRefSampleGeneExpIndexMap.get(refSample.Id);

                if(refSampleIndex == null)
                    continue;

                final double[] otherSampleTPMs = mRefSampleGeneExpression.getRow(refSampleIndex);

                double css = calcCosineSim(sampleTPMs, otherSampleTPMs);

                if(css < GENE_EXP_CSS_THRESHOLD)
                    continue;

                if(mConfig.WriteSimilarities)
                {
                    recordCssSimilarity(
                            topMatches, sample.Id, refSample.Id, css, EXPRESSION_PAIRWISE.toString(),
                            CSS_SIMILARITY_MAX_MATCHES, CSS_SIMILARITY_CUTOFF);
                }

                double cssWeight = pow(mCssExponent, -100 * (1 - css));

                int cancerSampleCount = mRefCancerSampleCounts.get(formCancerType(refSample));
                double weightedCss = css * cssWeight / sqrt(cancerSampleCount);

                totalWeightedCss += weightedCss;
            }

            cancerCssTotals.put(refCancerType, totalWeightedCss);
        }

        double totalCss = cancerCssTotals.values().stream().mapToDouble(x -> x).sum();

        for(Map.Entry<String,Double> entry : cancerCssTotals.entrySet())
        {
            double prob = totalCss > 0 ? entry.getValue() / totalCss : 0;
            cancerCssTotals.put(entry.getKey(), prob);
        }

        results.add(new SampleResult(
                sample.Id, GENE_EXP, CLASSIFIER, EXPRESSION_PAIRWISE.toString(), format("%.4g", totalCss), cancerCssTotals));

        similarities.addAll(topMatches);
    }

    private void addCancerCssResults(final SampleData sample, final double[] sampleGeneTPMs, final List<SampleResult> results)
    {
        int refCancerCount = mRefCancerTypeGeneExpression.Rows;

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
                    adjustRefCounts(mRefCancerTypeGeneExpression.getRow(i), sampleGeneTPMs, 1) : mRefCancerTypeGeneExpression.getRow(i);

            double css = calcCosineSim(sampleGeneTPMs, refPosFreqs);

            if(css < GENE_EXP_CSS_THRESHOLD)
                continue;

            double cssWeight = pow(mCssExponent, -100 * (1 - css));

            double weightedCss = css * cssWeight;
            cancerCssTotals.put(refCancerType, weightedCss);
        }

        double totalCss = cancerCssTotals.values().stream().mapToDouble(x -> x).sum();

        for(Map.Entry<String,Double> entry : cancerCssTotals.entrySet())
        {
            cancerCssTotals.put(entry.getKey(), entry.getValue() / totalCss);
        }

        results.add(new SampleResult(
                sample.Id, GENE_EXP, CLASSIFIER, EXPRESSION_COHORT.toString(), format("%.4g", totalCss), cancerCssTotals));
    }
}
