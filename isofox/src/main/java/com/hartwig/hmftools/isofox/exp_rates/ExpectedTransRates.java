package com.hartwig.hmftools.isofox.exp_rates;

import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.exp_rates.ExpectedRatesGenerator.formTranscriptDefinitions;
import static com.hartwig.hmftools.isofox.exp_rates.ExpectedRatesGenerator.writeExpectedRates;
import static com.hartwig.hmftools.common.sigs.DataUtils.RESIDUAL_PERC;
import static com.hartwig.hmftools.common.sigs.DataUtils.RESIDUAL_TOTAL;
import static com.hartwig.hmftools.common.sigs.DataUtils.calcResiduals;
import static com.hartwig.hmftools.common.sigs.DataUtils.calculateFittedCounts;
import static com.hartwig.hmftools.common.sigs.DataUtils.sumVector;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.gc.GcRatioCounts;
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
        double fitTotal = sumVector(fitAllocations);

        double[] residuals = calcResiduals(transComboCounts, fittedCounts, totalCounts);

        ISF_LOGGER.debug(String.format("gene(%s) totalFragments(%.0f) fitTotal(%.0f) residuals(%.0f perc=%.3f)",
                geneSummaryData.GeneNames, totalCounts, fitTotal, residuals[RESIDUAL_TOTAL], residuals[RESIDUAL_PERC]));

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
            writeCategoryCounts(mResultsWriter.getCategoryCountsWriter(), geneSummaryData.ChrId, mCurrentExpRatesData.Categories,
                    geneSummaryData.TransCategoryCounts, transComboCounts, fittedCounts, mConfig.ApplyGcBiasAdjust);
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
            double fragmentCount = tcData.fragmentCount();

            if(fragmentCount > 0)
            {
                int categoryId = mCurrentExpRatesData.getCategoryIndex(categoryKey);

                // for now if a category isn't found just log and then ignore the count in it
                if(categoryId < 0)
                {
                    ISF_LOGGER.trace("category({}) {} fragCount({}) skipped",
                            categoryKey, tcData.impliedType(), String.format("%.2f", fragmentCount));
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

    public static BufferedWriter createWriter(final IsofoxConfig config)
    {
        if(config.OutputDir.isEmpty())
            return null;

        try
        {
            final String outputFileName = config.formOutputFile("category_counts.csv");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("GenesId,Category,Count,FitCount");

            if(config.ApplyGcBiasAdjust)
            {
                GcRatioCounts tmp = new GcRatioCounts();

                for(Double gcRatio : tmp.getRatios())
                {
                    writer.write(String.format(",Gcr_%.2f", gcRatio));
                }
            }

            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write category counts data file: {}", e.toString());
            return null;
        }
    }

    private synchronized static void writeCategoryCounts(final BufferedWriter writer, final String genesId,
            final List<String> categories, final List<CategoryCountsData> categoryCountsData,
            final double[] counts, final double[] fittedCounts, boolean writeGcData)
    {
        try
        {
            final GcRatioCounts tmp = new GcRatioCounts();

            for(int i = 0; i < categories.size(); ++i)
            {
                double count = counts[i];
                final String category = categories.get(i);

                writer.write(String.format("%s,%s,%.0f,%.1f",
                        genesId, category, count, fittedCounts[i]));

                if(writeGcData)
                {
                    final CategoryCountsData catCounts = categoryCountsData.stream()
                            .filter(x -> x.combinedKey().equals(category)).findFirst().orElse(null);

                    final double[] gcCounts = catCounts != null ? catCounts.fragmentCountsByGcRatio() : tmp.getCounts();

                    for(int j = 0; j < gcCounts.length; ++j)
                    {
                        writer.write(String.format(",%.2f", gcCounts[j]));
                    }
                }

                writer.newLine();
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write category counts data file: {}", e.toString());
        }
    }

}
