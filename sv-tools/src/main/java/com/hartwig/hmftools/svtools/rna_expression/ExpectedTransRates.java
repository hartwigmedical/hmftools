package com.hartwig.hmftools.svtools.rna_expression;

import static com.hartwig.hmftools.sig_analyser.common.DataUtils.RESIDUAL_PERC;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.RESIDUAL_TOTAL;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.calcResiduals;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.calculateFittedCounts;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.sumVector;
import static com.hartwig.hmftools.svtools.rna_expression.ExpectedRatesData.ER_COL_CAT;
import static com.hartwig.hmftools.svtools.rna_expression.ExpectedRatesData.ER_COL_GENE_ID;
import static com.hartwig.hmftools.svtools.rna_expression.ExpectedRatesData.ER_COL_RATE;
import static com.hartwig.hmftools.svtools.rna_expression.ExpectedRatesData.ER_COL_TRANS_NAME;
import static com.hartwig.hmftools.svtools.rna_expression.ExpectedRatesGenerator.UNSPLICED_CAT_INDEX;
import static com.hartwig.hmftools.svtools.rna_expression.ExpectedRatesGenerator.formCategory;
import static com.hartwig.hmftools.svtools.rna_expression.FragmentMatchType.SHORT;
import static com.hartwig.hmftools.svtools.rna_expression.FragmentMatchType.SPLICED;
import static com.hartwig.hmftools.svtools.rna_expression.GeneMatchType.typeAsInt;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpConfig.RE_LOGGER;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class ExpectedTransRates
{
    private final RnaExpConfig mConfig;
    private final ResultsWriter mResultsWriter;

    private final Map<String, ExpectedRatesData> mGeneExpDataMap;

    private ExpectedRatesData mCurrentExpRatesData;

    public static final String UNSPLICED_ID = "UNSPLICED";

    public ExpectedTransRates(final RnaExpConfig config, final ResultsWriter resultsWriter)
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

    public static ExpectedTransRates from(final RnaExpConfig config)
    {
        return new ExpectedTransRates(config, null);
    }

    public boolean validData()
    {
        return mCurrentExpRatesData != null && mCurrentExpRatesData.validData();
    }

    public void loadGeneExpectedRatesData(final GeneReadData geneReadData)
    {
        mCurrentExpRatesData = mGeneExpDataMap.get(geneReadData.GeneData.GeneId);

        if(mCurrentExpRatesData == null)
        {
            RE_LOGGER.warn("gene({}) expected rates data not loaded", geneReadData.name());
        }
    }

    public void runTranscriptEstimation(
            final GeneReadData geneReadData, final List<TranscriptComboData> transCategoryCounts)
    {
        loadGeneExpectedRatesData(geneReadData);

        if(!validData())
        {
            RE_LOGGER.debug("gene({}) invalid expected rates or actuals data", geneReadData.name());
            return;
        }

        final double[] transComboCounts = generateTranscriptCounts(geneReadData, transCategoryCounts);

        double totalCounts = sumVector(transComboCounts);

        if(totalCounts == 0)
            return;

        // add in counts for the unspliced category
        int unsplicedCount = geneReadData.getCounts()[typeAsInt(GeneMatchType.UNSPLICED)];
        transComboCounts[UNSPLICED_CAT_INDEX] = unsplicedCount;

        final List<String> transcriptNames = mCurrentExpRatesData.TranscriptIds;

        final double[] fitAllocations = ExpectationMaxFit.performFit(transComboCounts, mCurrentExpRatesData.getTranscriptDefinitions());
        final double[] fittedCounts = calculateFittedCounts(mCurrentExpRatesData.getTranscriptDefinitions(), fitAllocations);

        double[] residuals = calcResiduals(transComboCounts, fittedCounts, totalCounts);

        RE_LOGGER.debug(String.format("gene(%s) totalFragments(%.0f) residuals(%.0f perc=%.3f)",
                geneReadData.name(), totalCounts, residuals[RESIDUAL_TOTAL], residuals[RESIDUAL_PERC]));

        geneReadData.setFitResiduals(residuals[RESIDUAL_TOTAL]);

        Map<String,Double> transAllocations = geneReadData.getTranscriptAllocations();

        for(int transIndex = 0; transIndex < transcriptNames.size(); ++transIndex)
        {
            double transAllocation = fitAllocations[transIndex];
            final String transName = transcriptNames.get(transIndex);

            if(transAllocation > 0)
            {
                RE_LOGGER.debug("transcript({}) allocated count({})", transName, String.format("%.2f", transAllocation));
            }

            transAllocations.put(transName, transAllocation);
        }

        if(mConfig.WriteTransComboData)
        {
            mResultsWriter.writeTransComboCounts(
                    geneReadData, mCurrentExpRatesData.Categories, transComboCounts, fittedCounts);
        }
    }

    public double[] generateTranscriptCounts(final GeneReadData geneReadData, final List<TranscriptComboData> transComboData)
    {
        double[] categoryCounts = new double[mCurrentExpRatesData.Categories.size()];

        int skippedComboCounts = 0;

        for(TranscriptComboData tcData : transComboData)
        {
            final String transKey = !tcData.getTranscriptIds().isEmpty() ? tcData.getTranscriptsKey() : UNSPLICED_ID;

            int shortCount = tcData.getShortCount();
            int splicedCount = tcData.getSplicedCount();

            if(shortCount > 0)
            {
                final String categoryStr = formCategory(transKey, SHORT);
                int categoryId = mCurrentExpRatesData.getCategoryIndex(categoryStr);

                // for now if a category isn't found just log and then ignore the count in it
                if(categoryId < 0)
                {
                    RE_LOGGER.debug("category({}) skipped with count({})", categoryStr, shortCount);
                    skippedComboCounts += shortCount;
                }
                else
                {
                    categoryCounts[categoryId] = shortCount; //  * mConfig.UnsplicedWeight;
                }
            }

            if(splicedCount > 0)
            {
                final String categoryStr = formCategory(transKey, SPLICED);
                int categoryId = mCurrentExpRatesData.getCategoryIndex(categoryStr);

                if(categoryId < 0)
                {
                    RE_LOGGER.debug("category({}) skipped with count({})", categoryStr, splicedCount);
                    skippedComboCounts += splicedCount;
                }
                else
                {
                    categoryCounts[categoryId] = splicedCount;
                }
            }
        }

        if(skippedComboCounts > 0)
        {
            double totalCounts = sumVector(categoryCounts) + skippedComboCounts;

            RE_LOGGER.debug(String.format("gene(%s) skippedCounts(%d perc=%.3f of total=%.0f)",
                    geneReadData.GeneData.GeneName, skippedComboCounts, skippedComboCounts/totalCounts, totalCounts));
        }

        return categoryCounts;
    }

    public static void loadExpRatesFile(final String filename, final Map<String, ExpectedRatesData> geneExpDataMap)
    {
        if (!Files.exists(Paths.get(filename)))
        {
            RE_LOGGER.warn("invalid gene ID file({})", filename);
            return;
        }

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            // skip field names
            String line = fileReader.readLine();

            if (line == null)
            {
                RE_LOGGER.error("Empty patient sample IDs file({})", filename);
                return;
            }

            List<String[]> geneRatesData = Lists.newArrayList();
            ExpectedRatesData expExpData = null;

            while ((line = fileReader.readLine()) != null)
            {
                String[] items = line.split(",");

                if (items.length != ER_COL_RATE + 1)
                {
                    RE_LOGGER.error("invalid exp rates data: {}", line);
                    return;
                }

                String geneId = items[ER_COL_GENE_ID];
                String transName = items[ER_COL_TRANS_NAME];
                String categoryStr = items[ER_COL_CAT];

                if(expExpData == null)
                {
                    expExpData = new ExpectedRatesData(geneId);
                }
                else if(!expExpData.GeneId.equals(geneId))
                {
                    expExpData.buildDefinitionsFromFileData(geneRatesData);
                    geneExpDataMap.put(expExpData.GeneId, expExpData);

                    // start the new one
                    geneRatesData.clear();
                    expExpData = new ExpectedRatesData(geneId);
                }

                expExpData.addCategory(categoryStr);

                if(!expExpData.TranscriptIds.contains(transName))
                    expExpData.TranscriptIds.add(transName);

                geneRatesData.add(items);
            }

            RE_LOGGER.info("loaded {} gene expected trans exp rates from file({})", geneExpDataMap.size(), filename);
        }
        catch (IOException e)
        {
            RE_LOGGER.warn("failed to load expected rates file({}): {}", filename, e.toString());
        }
    }


}
