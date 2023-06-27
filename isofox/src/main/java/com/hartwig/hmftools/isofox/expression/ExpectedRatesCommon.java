package com.hartwig.hmftools.isofox.expression;

import static com.hartwig.hmftools.common.sigs.SigUtils.convertToPercentages;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

public final class ExpectedRatesCommon
{
    public static final int FL_LENGTH = 0;
    public static final int FL_FREQUENCY = 1;
    public static final String EXP_COUNT_LENGTH_HEADER = "Length_";

    public static Map<String,List<CategoryCountsData>> createTransComboDataMap(final List<CategoryCountsData> categoryCountsData)
    {
        final Map<String,List<CategoryCountsData>> transGeneCountsMap = Maps.newHashMap();

        for(final CategoryCountsData tcData : categoryCountsData)
        {
            final Set<String> geneTransNames = Sets.newHashSet();
            tcData.unsplicedGeneIds().forEach(x -> geneTransNames.add(x));
            tcData.transcriptIds().forEach(x -> geneTransNames.add(String.valueOf(x)));

            for(String geneTransName : geneTransNames)
            {
                List<CategoryCountsData> countsList = transGeneCountsMap.get(geneTransName);
                if(countsList == null)
                {
                    transGeneCountsMap.put(geneTransName, Lists.newArrayList(tcData));
                }
                else
                {
                    countsList.add(tcData);
                }
            }
        }

        return transGeneCountsMap;
    }

    public static void formTranscriptDefinitions(final List<CategoryCountsData> categoryCountsData, ExpectedRatesData expRatesData)
    {
        // convert fragment counts in each category per transcript into the equivalent of a signature per transcript
        collectCategories(categoryCountsData, expRatesData);

        final Map<String,List<CategoryCountsData>> transGeneCountsMap = createTransComboDataMap(categoryCountsData);

        int categoryCount = expRatesData.Categories.size();

        transGeneCountsMap.keySet().forEach(x -> expRatesData.TranscriptIds.add(x));

        expRatesData.initialiseTranscriptDefinitions();

        for(int transIndex = 0; transIndex < expRatesData.TranscriptIds.size(); ++transIndex)
        {
            final String transId = expRatesData.TranscriptIds.get(transIndex);

            double[] categoryCounts = new double[categoryCount];

            final List<CategoryCountsData> transCounts = transGeneCountsMap.get(transId);

            for(CategoryCountsData tcData : transCounts)
            {
                final String transKey = tcData.combinedKey();
                double fragmentCount = tcData.fragmentCount();

                if(fragmentCount > 0)
                {
                    int categoryId = expRatesData.getCategoryIndex(transKey);

                    if(categoryId < 0)
                    {
                        ISF_LOGGER.error("invalid category index from transKey({})", transKey);
                        return;
                    }

                    categoryCounts[categoryId] = fragmentCount;
                }
            }

            // convert counts to ratios and add against this transcript definition
            convertToPercentages(categoryCounts);
            expRatesData.getTranscriptDefinitions().setCol(transIndex, categoryCounts);
        }
    }

    private static void collectCategories(
            final List<CategoryCountsData> categoryCountsData, ExpectedRatesData expRatesData)
    {
        for(final CategoryCountsData tcData : categoryCountsData)
        {
            final String transKey = tcData.combinedKey();

            if(tcData.fragmentCount() > 0 || tcData.transcriptIds().isEmpty()) // force inclusion of unspliced gene categories
            {
                expRatesData.addCategory(transKey);
            }
        }
    }
}
