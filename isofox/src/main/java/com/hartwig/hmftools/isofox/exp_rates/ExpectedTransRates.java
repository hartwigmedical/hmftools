package com.hartwig.hmftools.isofox.exp_rates;

import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.exp_rates.ExpectedRatesGenerator.formTranscriptDefinitions;
import static com.hartwig.hmftools.isofox.exp_rates.ExpectedRatesGenerator.writeExpectedRates;
import static com.hartwig.hmftools.common.sigs.DataUtils.RESIDUAL_PERC;
import static com.hartwig.hmftools.common.sigs.DataUtils.RESIDUAL_TOTAL;
import static com.hartwig.hmftools.common.sigs.DataUtils.calcResiduals;
import static com.hartwig.hmftools.common.sigs.DataUtils.calculateFittedCounts;
import static com.hartwig.hmftools.common.sigs.DataUtils.sumVector;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.results.GeneResult;
import com.hartwig.hmftools.isofox.results.ResultsWriter;

public class ExpectedTransRates
{
    private final IsofoxConfig mConfig;
    private final ResultsWriter mResultsWriter;
    private final ExpectedCountsCache mCache;

    private ExpectedRatesData mCurrentExpRatesData;

    public ExpectedTransRates(final IsofoxConfig config, final ExpectedCountsCache cache, final ResultsWriter resultsWriter)
    {
        mConfig = config;
        mResultsWriter = resultsWriter;

        mCache = cache;
        mCurrentExpRatesData = null;
    }

    public static ExpectedTransRates from(final IsofoxConfig config)
    {
        return new ExpectedTransRates(config, null,null);
    }

    public boolean validData()
    {
        return mCurrentExpRatesData != null && mCurrentExpRatesData.validData();
    }

    private void applyFragmentLengthDistributionToExpectedCounts(final Map<String,List<CategoryCountsData>> geneSetCountsData)
    {
        final List<int[]> fragmentLengthData = mConfig.FragmentLengthData;

        for(List<CategoryCountsData> categoryCountsData : geneSetCountsData.values())
        {
            categoryCountsData.forEach(x -> x.applyFrequencies(fragmentLengthData));
        }
    }

    public void runTranscriptEstimation(final GeneCollectionSummaryData geneSummaryData, final ExpectedRatesData expRatesData)
    {
        if(expRatesData == null)
        {
            loadGeneExpectedRatesData(geneSummaryData.ChrId, geneSummaryData.GeneIds);
        }
        else
        {
            mCurrentExpRatesData = expRatesData;
        }

        if(!validData())
        {
            ISF_LOGGER.debug("gene({}) invalid expected rates or actuals data", geneSummaryData.GeneNames);
            return;
        }

        final double[] transComboCounts = generateReadCounts(geneSummaryData);

        double totalCounts = sumVector(transComboCounts);

        if(totalCounts == 0)
            return;

        final List<String> transcriptNames = mCurrentExpRatesData.TranscriptIds;

        final double[] fitAllocations = ExpectationMaxFit.performFit(transComboCounts, mCurrentExpRatesData.getTranscriptDefinitions());
        final double[] fittedCounts = calculateFittedCounts(mCurrentExpRatesData.getTranscriptDefinitions(), fitAllocations);

        double[] residuals = calcResiduals(transComboCounts, fittedCounts, totalCounts);

        ISF_LOGGER.debug(String.format("gene(%s) totalFragments(%.0f) residuals(%.0f perc=%.3f)",
                geneSummaryData.GeneNames, totalCounts, residuals[RESIDUAL_TOTAL], residuals[RESIDUAL_PERC]));

        geneSummaryData.setFitResiduals(residuals[RESIDUAL_TOTAL]);

        final Map<String,Double> transAllocations = geneSummaryData.getFitAllocations();

        for(int transIndex = 0; transIndex < transcriptNames.size(); ++transIndex)
        {
            double transAllocation = fitAllocations[transIndex];
            final String transName = transcriptNames.get(transIndex);

            if(transAllocation > 0)
            {
                ISF_LOGGER.trace("transcript({}) allocated count({})", transName, String.format("%.2f", transAllocation));
            }

            transAllocations.put(transName, transAllocation);
        }

        if(mConfig.WriteTransComboData)
        {
            mResultsWriter.writeTransComboCounts(
                    geneSummaryData.ChrId, mCurrentExpRatesData.Categories, transComboCounts, fittedCounts);
        }
    }

    public static double calcTotalTranscriptExpression(final List<GeneCollectionSummaryData> geneSummaryData)
    {
        double totalFragsPerKb = 0;

        for(final GeneCollectionSummaryData summaryData : geneSummaryData)
        {
            totalFragsPerKb += summaryData.TranscriptResults.stream().mapToDouble(x -> x.fragmentsPerKb()).sum();
        }

        return totalFragsPerKb;
    }

    public static void setTranscriptsPerMillion(final List<GeneCollectionSummaryData> geneSummaryData, double tpmFactor)
    {
        for(final GeneCollectionSummaryData summaryData : geneSummaryData)
        {
            summaryData.TranscriptResults.forEach(x -> x.setTPM(x.fragmentsPerKb()/tpmFactor));
        }
    }

    private void loadGeneExpectedRatesData(final String chrId, final List<String> geneIds)
    {
        mCurrentExpRatesData = null;

        if(mCache.hasExpectedRatesCached())
        {
            mCurrentExpRatesData = mCache.getGeneExpectedRatesData(chrId);
            return;
        }

        final Map<String,List<CategoryCountsData>> geneSetCountsData = mCache.getGeneExpectedRatesData(chrId, geneIds);

        if(geneSetCountsData == null)
        {
            ISF_LOGGER.warn("genes({}: {}) expected counts data not loaded", chrId, appendStrList(geneIds, ';'));
            return;
        }

        mCurrentExpRatesData = new ExpectedRatesData(chrId);

        // apply observed fragment length distribution to the generated counts
        if(mConfig.ApplyFragmentLengthAdjust)
            applyFragmentLengthDistributionToExpectedCounts(geneSetCountsData);

        formTranscriptDefinitions(geneSetCountsData, mCurrentExpRatesData);

        if(mConfig.WriteExpectedRates)
        {
            writeExpectedRates(mResultsWriter.getExpRatesWriter(), mCurrentExpRatesData);
        }
    }

    private double[] generateReadCounts(final GeneCollectionSummaryData geneSummaryData)
    {
        double[] categoryCounts = new double[mCurrentExpRatesData.Categories.size()];

        int skippedComboCounts = 0;

        for(CategoryCountsData tcData : geneSummaryData.TransCategoryCounts)
        {
            final String categoryKey = tcData.combinedKey();
            long fragmentCount = tcData.fragmentCount();

            if(fragmentCount > 0)
            {
                int categoryId = mCurrentExpRatesData.getCategoryIndex(categoryKey);

                // for now if a category isn't found just log and then ignore the count in it
                if(categoryId < 0)
                {
                    ISF_LOGGER.trace("category({}) {} fragCount({}) skipped", categoryKey, tcData.impliedType(), fragmentCount);
                    skippedComboCounts += fragmentCount;
                }
                else
                {
                    categoryCounts[categoryId] = fragmentCount;
                }
            }
        }

        if(skippedComboCounts > 0)
        {
            double totalCounts = sumVector(categoryCounts) + skippedComboCounts;

            ISF_LOGGER.debug(String.format("gene(%s) skippedCounts(%d perc=%.3f of total=%.0f)",
                    geneSummaryData.GeneNames, skippedComboCounts, skippedComboCounts/totalCounts, totalCounts));
        }

        return categoryCounts;
    }

}
