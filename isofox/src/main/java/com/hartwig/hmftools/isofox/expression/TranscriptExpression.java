package com.hartwig.hmftools.isofox.expression;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.sigs.SigUtils.calcResiduals;
import static com.hartwig.hmftools.common.sigs.SigUtils.calculateFittedCounts;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;
import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.MAX_GENE_PERC_CONTRIBUTION;
import static com.hartwig.hmftools.isofox.expression.CategoryCountsData.hasGeneIdentifier;
import static com.hartwig.hmftools.isofox.expression.ExpectedRatesGenerator.formTranscriptDefinitions;
import static com.hartwig.hmftools.isofox.expression.ExpectedRatesGenerator.writeExpectedRates;

import static org.apache.logging.log4j.Level.DEBUG;
import static org.apache.logging.log4j.Level.WARN;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sigs.ExpectationMaxFit;
import com.hartwig.hmftools.common.sigs.SigResiduals;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.adjusts.FragmentSize;
import com.hartwig.hmftools.isofox.adjusts.GcRatioCounts;
import com.hartwig.hmftools.isofox.results.GeneResult;
import com.hartwig.hmftools.isofox.results.ResultsWriter;
import com.hartwig.hmftools.isofox.results.TranscriptResult;

import org.apache.logging.log4j.Level;

public class TranscriptExpression
{
    private final IsofoxConfig mConfig;
    private final ResultsWriter mResultsWriter;
    private final ExpectedCountsCache mCache;

    private final Map<String,ExpectedRatesData> mExpectedRatesDataMap; // cached computed results for this chromosome
    private ExpectedRatesData mCurrentExpRatesData;

    private final List<Double> mFragmentFrequencyRates;

    public TranscriptExpression(final IsofoxConfig config, final ExpectedCountsCache cache, final ResultsWriter resultsWriter)
    {
        mConfig = config;
        mResultsWriter = resultsWriter;

        mCache = cache;
        mExpectedRatesDataMap = Maps.newHashMap();
        mCurrentExpRatesData = null;

        // convert fragment distribution counts to rates
        double totalFrequencyFrags = mConfig.FragmentSizeData.stream().mapToDouble(x -> x.Frequency).sum();
        mFragmentFrequencyRates = Lists.newArrayList();

        if(totalFrequencyFrags > 0)
        {
            for(final FragmentSize fragFrequency : mConfig.FragmentSizeData)
            {
                mFragmentFrequencyRates.add(fragFrequency.Frequency / totalFrequencyFrags);
            }
        }
    }

    public static TranscriptExpression from(final IsofoxConfig config)
    {
        return new TranscriptExpression(config, null,null);
    }

    public boolean validData()
    {
        return mCurrentExpRatesData != null && mCurrentExpRatesData.validData();
    }

    private void applyFragmentLengthDistributionToExpectedCounts(final List<CategoryCountsData> geneSetCountsData)
    {
        geneSetCountsData.forEach(x -> x.applyFrequencies(mFragmentFrequencyRates));
    }

    public void runTranscriptEstimation(
            final Map<Integer,String> transIdMap, final GeneCollectionSummary geneSummaryData,
            final ExpectedRatesData expRatesData, boolean checkCached)
    {
        if(expRatesData == null)
        {
            if(checkCached)
                mCurrentExpRatesData = mExpectedRatesDataMap.get(geneSummaryData.ChrId);
            else
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

        if(!checkCached && mConfig.applyGcBiasAdjust()) // cache the generated data since it will be used again in GC adjustment calcs
            mExpectedRatesDataMap.put(geneSummaryData.ChrId, mCurrentExpRatesData);

        final double[] transComboCounts = generateReadCounts(geneSummaryData);

        double totalCounts = sumVector(transComboCounts);

        if(totalCounts == 0)
            return;

        final List<String> transcriptIds = mCurrentExpRatesData.TranscriptIds;

        final double[] fitAllocations = ExpectationMaxFit.performFit(transComboCounts, mCurrentExpRatesData.getTranscriptDefinitions());
        final double[] fittedCounts = calculateFittedCounts(mCurrentExpRatesData.getTranscriptDefinitions(), fitAllocations);
        double fitTotal = sumVector(fitAllocations);

        SigResiduals residuals = calcResiduals(transComboCounts, fittedCounts, totalCounts);

        ISF_LOGGER.debug(String.format("gene(%s) totalFragments(%.0f) fitTotal(%.0f) residuals(%.0f perc=%.3f)",
                geneSummaryData.GeneNames, totalCounts, fitTotal, residuals.Total, residuals.Percent));

        geneSummaryData.setFitResiduals(residuals.Total);

        final Map<String,Double> transAllocations = geneSummaryData.getFitAllocations();

        for(int transIndex = 0; transIndex < transcriptIds.size(); ++transIndex)
        {
            final String transGeneId = transcriptIds.get(transIndex);
            final String transName = hasGeneIdentifier(transGeneId) ? transGeneId : transIdMap.get(Integer.parseInt(transGeneId));

            double transAllocation = fitAllocations[transIndex];

            if(transAllocation > 0)
            {
                ISF_LOGGER.trace("transcript({}) allocated count({})", transName, String.format("%.2f", transAllocation));
            }

            transAllocations.put(transName, transAllocation);
        }

        if(mConfig.WriteTransComboData)
        {
            writeCategoryCounts(mResultsWriter.getCategoryCountsWriter(), geneSummaryData.ChrId, mCurrentExpRatesData.Categories,
                    geneSummaryData.TransCategoryCounts, transComboCounts, fittedCounts, mConfig.applyGcBiasAdjust());
        }
    }

    private static final int RAW_TPM = 0;
    private static final int ADJUSTED_TPM = 1;
    private static final double TPM_MILLION = 1000000;

    public static double[] calcTpmFactors(final List<GeneCollectionSummary> geneSummaryData, final List<String> enrichedGeneIds)
    {
        // exclude enriched genes and cap the contribution of any one gene to 1%
        double[] results = {0, 0};

        List<Double> fragsPerKbSet = Lists.newArrayListWithExpectedSize(200000);

        double maxFragPerKb = 0;
        double fragsPerKbTotal = 0;
        double rawFragsPerKbTotal = 0;

        for(final GeneCollectionSummary summaryData : geneSummaryData)
        {
            for(final TranscriptResult transResult : summaryData.TranscriptResults)
            {
                double fragsPerKb = transResult.fragmentsPerKb();

                rawFragsPerKbTotal += fragsPerKb;

                if(enrichedGeneIds.contains(transResult.Trans.GeneId))
                    continue;

                if(fragsPerKb < 0.5)
                {
                    fragsPerKbTotal += fragsPerKb;
                }
                else
                {
                    fragsPerKbSet.add(fragsPerKb);
                    fragsPerKbTotal += fragsPerKb;
                    maxFragPerKb = max(maxFragPerKb, fragsPerKb);
                }
            }
        }

        results[RAW_TPM] = rawFragsPerKbTotal / TPM_MILLION;

        if(maxFragPerKb <= MAX_GENE_PERC_CONTRIBUTION * fragsPerKbTotal)
        {
            results[ADJUSTED_TPM] = fragsPerKbTotal / TPM_MILLION;
        }
        else
        {
            int iterations = 0;
            int maxIterations = 10;

            double maxFragsPerKbTotal = fragsPerKbTotal;
            double nextFragsPerKbTotal = fragsPerKbTotal * 0.5;
            double minFragsPerKbTotal = 0;

            // frag total is initially too high so halve it (the new min) and set the max to this value
            while (iterations < maxIterations)
            {
                double maxAllowedContrib = nextFragsPerKbTotal * MAX_GENE_PERC_CONTRIBUTION;
                double maxFragValue = 0;
                fragsPerKbTotal = 0;
                for(Double fragsPerKb : fragsPerKbSet)
                {
                    double fragValue = min(fragsPerKb, maxAllowedContrib);
                    maxFragValue = max(maxFragValue, fragValue);
                    fragsPerKbTotal += fragValue;
                }

                double maxGeneContribPerc = maxFragValue / fragsPerKbTotal;

                if(abs(maxGeneContribPerc - MAX_GENE_PERC_CONTRIBUTION) < 0.001)
                    break;

                if (maxGeneContribPerc > MAX_GENE_PERC_CONTRIBUTION)
                {
                    maxFragsPerKbTotal = fragsPerKbTotal;
                    nextFragsPerKbTotal = max(fragsPerKbTotal * 0.5, minFragsPerKbTotal);
                }
                else
                {
                    // frag total was too low
                    minFragsPerKbTotal = fragsPerKbTotal;
                    nextFragsPerKbTotal = (fragsPerKbTotal + maxFragsPerKbTotal) * 0.5;
                }

                ++iterations;
            }
        }

        results[ADJUSTED_TPM] = fragsPerKbTotal / TPM_MILLION;
        return results;
    }

    public static void setTranscriptsPerMillion(
            final List<GeneCollectionSummary> geneSummaryData, double[] tmpFactors)
    {
        for(final GeneCollectionSummary summaryData : geneSummaryData)
        {
            Map<String,double[]> geneTPMs = Maps.newHashMap();

            for(TranscriptResult transResult : summaryData.TranscriptResults)
            {
                double fragsPerKb = transResult.fragmentsPerKb();

                double rawTpm = tmpFactors[RAW_TPM] > 0 ? fragsPerKb/tmpFactors[RAW_TPM] : 0;
                double adjustedTpm = tmpFactors[ADJUSTED_TPM] > 0 ? fragsPerKb/tmpFactors[ADJUSTED_TPM] : 0;
                transResult.setTPM(rawTpm, adjustedTpm);

                double[] geneTpm = geneTPMs.get(transResult.Trans.GeneId);
                if(geneTpm == null)
                {
                    geneTPMs.put(transResult.Trans.GeneId, new double[] { rawTpm, adjustedTpm });
                }
                else
                {
                    geneTpm[RAW_TPM] += rawTpm;
                    geneTpm[ADJUSTED_TPM] += adjustedTpm;
                }
            }

            for(GeneResult geneResult : summaryData.GeneResults)
            {
                final double[] geneTpm = geneTPMs.get(geneResult.GeneData.GeneId);
                geneResult.setTPM(geneTpm[RAW_TPM], geneTpm[ADJUSTED_TPM]);
            }
        }
    }

    private void loadGeneExpectedRatesData(final String chrId, final List<String> geneIds)
    {
        mCurrentExpRatesData = null;

        final List<CategoryCountsData> geneSetCountsData = mCache.getGeneExpectedRatesData(chrId, geneIds);

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

    private double[] generateReadCounts(final GeneCollectionSummary geneSummaryData)
    {
        double[] categoryCounts = new double[mCurrentExpRatesData.Categories.size()];

        int skippedComboCounts = 0;
        double totalCounts = geneSummaryData.TransCategoryCounts.stream().mapToDouble(x -> x.fragmentCount()).sum();

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
                    if(fragmentCount/totalCounts > 0.20 && fragmentCount > 50)
                    {
                        ISF_LOGGER.debug("category({}) {} fragCount({}) skipped",
                                categoryKey, tcData.impliedType(), String.format("%.2f", fragmentCount));
                    }
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
            double skippedPerc = skippedComboCounts/totalCounts;

            ISF_LOGGER.debug(String.format("gene(%d:%s) categories(act=%d exp=%d trans+genes=%d) skippedCounts(%d perc=%.3f of total=%.0f)",
                    geneSummaryData.GeneIds.size(), geneSummaryData.GeneNames,
                    geneSummaryData.TransCategoryCounts.size(), mCurrentExpRatesData.Categories.size(),
                    mCurrentExpRatesData.TranscriptIds.size(), skippedComboCounts, skippedPerc, totalCounts));
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

            if(config.applyGcBiasAdjust())
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
