package com.hartwig.hmftools.isofox.exp_rates;

import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.FragmentMatchType.SHORT;
import static com.hartwig.hmftools.isofox.common.FragmentMatchType.SPLICED;
import static com.hartwig.hmftools.isofox.exp_rates.ExpectedRatesData.ER_COL_GENE_SET_ID;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.RESIDUAL_PERC;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.RESIDUAL_TOTAL;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.calcResiduals;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.calculateFittedCounts;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.sumVector;
import static com.hartwig.hmftools.isofox.exp_rates.ExpectedRatesData.ER_COL_CAT;
import static com.hartwig.hmftools.isofox.exp_rates.ExpectedRatesData.ER_COL_RATE;
import static com.hartwig.hmftools.isofox.exp_rates.ExpectedRatesData.ER_COL_TRANS_NAME;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.results.ResultsWriter;

public class ExpectedTransRates
{
    private final IsofoxConfig mConfig;
    private final ResultsWriter mResultsWriter;

    private final Map<String, ExpectedRatesData> mGeneExpDataMap;

    private ExpectedRatesData mCurrentExpRatesData;

    public ExpectedTransRates(final IsofoxConfig config, final ResultsWriter resultsWriter)
    {
        mConfig = config;
        mResultsWriter = resultsWriter;

        mGeneExpDataMap = Maps.newHashMap();
        mCurrentExpRatesData = null;

        if(mConfig.ExpRatesFile != null && Files.exists(Paths.get(mConfig.ExpRatesFile)))
        {
            loadExpRatesFile(mConfig.ExpRatesFile, mGeneExpDataMap);
        }
    }

    public static ExpectedTransRates from(final IsofoxConfig config)
    {
        return new ExpectedTransRates(config, null);
    }

    public boolean validData()
    {
        return mCurrentExpRatesData != null && mCurrentExpRatesData.validData();
    }

    public void loadGeneExpectedRatesData(final GeneCollection genes)
    {
        mCurrentExpRatesData = mGeneExpDataMap.get(genes.chrId());

        if(mCurrentExpRatesData == null)
        {
            ISF_LOGGER.warn("gene({}) expected rates data not loaded", genes.geneNames());
        }
    }

    public void runTranscriptEstimation(
            final GeneCollection genes, final List<TranscriptComboData> transCategoryCounts, final ExpectedRatesData expRatesData)
    {
        if(expRatesData == null)
        {
            loadGeneExpectedRatesData(genes);
        }
        else
        {
            mCurrentExpRatesData = expRatesData;
        }

        if(!validData())
        {
            ISF_LOGGER.debug("gene({}) invalid expected rates or actuals data", genes.geneNames());
            return;
        }

        final double[] transComboCounts = generateReadCounts(genes, transCategoryCounts);

        double totalCounts = sumVector(transComboCounts);

        if(totalCounts == 0)
            return;

        // add in counts for the unspliced category
        // int unsplicedCount = geneReadData.getCounts()[typeAsInt(GeneMatchType.UNSPLICED)];
        // transComboCounts[UNSPLICED_CAT_INDEX] = unsplicedCount;

        final List<String> transcriptNames = mCurrentExpRatesData.TranscriptIds;

        final double[] fitAllocations = ExpectationMaxFit.performFit(transComboCounts, mCurrentExpRatesData.getTranscriptDefinitions());
        final double[] fittedCounts = calculateFittedCounts(mCurrentExpRatesData.getTranscriptDefinitions(), fitAllocations);

        double[] residuals = calcResiduals(transComboCounts, fittedCounts, totalCounts);

        ISF_LOGGER.debug(String.format("gene(%s) totalFragments(%.0f) residuals(%.0f perc=%.3f)",
                genes.geneNames(), totalCounts, residuals[RESIDUAL_TOTAL], residuals[RESIDUAL_PERC]));

        genes.setFitResiduals(residuals[RESIDUAL_TOTAL]);

        Map<String,Double> transAllocations = genes.getTranscriptAllocations();

        for(int transIndex = 0; transIndex < transcriptNames.size(); ++transIndex)
        {
            double transAllocation = fitAllocations[transIndex];
            final String transName = transcriptNames.get(transIndex);

            if(transAllocation > 0)
            {
                ISF_LOGGER.debug("transcript({}) allocated count({})", transName, String.format("%.2f", transAllocation));
            }

            transAllocations.put(transName, transAllocation);
        }

        if(mConfig.WriteTransComboData)
        {
            mResultsWriter.writeTransComboCounts(
                    genes.chrId(), mCurrentExpRatesData.Categories, transComboCounts, fittedCounts);
        }
    }

    private double[] generateReadCounts(final GeneCollection genes, final List<TranscriptComboData> transComboData)
    {
        double[] categoryCounts = new double[mCurrentExpRatesData.Categories.size()];

        int skippedComboCounts = 0;

        for(TranscriptComboData tcData : transComboData)
        {
            final String categoryKey = tcData.combinedKey();
            int fragmentCount = tcData.fragmentCount();

            if(fragmentCount > 0)
            {
                int categoryId = mCurrentExpRatesData.getCategoryIndex(categoryKey);

                // for now if a category isn't found just log and then ignore the count in it
                if(categoryId < 0)
                {
                    ISF_LOGGER.debug("category({}) {} fragCount({}) skipped", categoryKey, tcData.impliedType(), fragmentCount);
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
                    genes.geneNames(), skippedComboCounts, skippedComboCounts/totalCounts, totalCounts));
        }

        return categoryCounts;
    }

    public static void loadExpRatesFile(final String filename, final Map<String, ExpectedRatesData> geneExpDataMap)
    {
        if (!Files.exists(Paths.get(filename)))
        {
            ISF_LOGGER.warn("invalid gene ID file({})", filename);
            return;
        }

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            // skip field names
            String line = fileReader.readLine();

            if (line == null)
            {
                ISF_LOGGER.error("Empty patient sample IDs file({})", filename);
                return;
            }

            List<String[]> geneRatesData = Lists.newArrayList();
            ExpectedRatesData expExpData = null;

            while ((line = fileReader.readLine()) != null)
            {
                String[] items = line.split(",");

                if (items.length != ER_COL_RATE + 1)
                {
                    ISF_LOGGER.error("invalid exp rates data: {}", line);
                    return;
                }

                String geneCollectionId = items[ER_COL_GENE_SET_ID];
                String transName = items[ER_COL_TRANS_NAME];
                String categoryStr = items[ER_COL_CAT];

                if(expExpData == null)
                {
                    expExpData = new ExpectedRatesData(geneCollectionId);
                }
                else if(!expExpData.Id.equals(geneCollectionId))
                {
                    expExpData.buildDefinitionsFromFileData(geneRatesData);
                    geneExpDataMap.put(expExpData.Id, expExpData);

                    // start the new one
                    geneRatesData.clear();
                    expExpData = new ExpectedRatesData(geneCollectionId);
                }

                expExpData.addCategory(categoryStr);

                if(!expExpData.TranscriptIds.contains(transName))
                    expExpData.TranscriptIds.add(transName);

                geneRatesData.add(items);
            }

            ISF_LOGGER.info("loaded {} gene expected trans exp rates from file({})", geneExpDataMap.size(), filename);
        }
        catch (IOException e)
        {
            ISF_LOGGER.warn("failed to load expected rates file({}): {}", filename, e.toString());
        }
    }


}
