package com.hartwig.hmftools.isofox.exp_rates;

import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.GeneCollection;

public class ExpectedCountsCache
{
    private final IsofoxConfig mConfig;
    private final Map<String,Map<String,List<CategoryCountsData>>> mGeneSetCategoryDataMap;
    private final boolean mValidData;

    public ExpectedCountsCache(final IsofoxConfig config)
    {
        mConfig = config;
        mGeneSetCategoryDataMap = Maps.newHashMap();

        if(config.ExpCountsFile != null && Files.exists(Paths.get(mConfig.ExpCountsFile)))
        {
            mValidData = loadExpCountsFile();
        }
        else
        {
            mValidData = false;
        }
    }

    public Map<String,List<CategoryCountsData>> getGeneExpectedRatesData(final GeneCollection genes)
    {
        Map<String, List<CategoryCountsData>> geneSetCountsData = mGeneSetCategoryDataMap.get(genes.chrId());

        if (geneSetCountsData == null || !geneSetCountsDataMatches(genes, geneSetCountsData.keySet()))
        {
            geneSetCountsData = findGeneSetCountsData(genes);
        }

        return geneSetCountsData;
    }

    private boolean geneSetCountsDataMatches(final GeneCollection geneCollection, final Set<String> geneTransSet)
    {
        // confirm that the genes in the collection match
        return !geneCollection.genes().stream().anyMatch(x -> !geneTransSet.contains(x.GeneData.GeneId));
    }

    private final Map<String,List<CategoryCountsData>> findGeneSetCountsData(final GeneCollection geneCollection)
    {
        for(final Map<String,List<CategoryCountsData>> geneCounts : mGeneSetCategoryDataMap.values())
        {
            if(geneSetCountsDataMatches(geneCollection, geneCounts.keySet()))
                return geneCounts;
        }

        return null;
    }

    // GeneCollectionId,TransId,Category,Counts for each fragment length
    public static final int ER_COL_GENE_SET_ID = 0;
    public static final int ER_COL_TRANS_GENE_NAME = 1;
    public static final int ER_COL_CAT = 2;

    private boolean loadExpCountsFile()
    {
        if (!Files.exists(Paths.get(mConfig.ExpCountsFile)))
        {
            ISF_LOGGER.warn("invalid gene ID file({})", mConfig.ExpCountsFile);
            return false;
        }

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(mConfig.ExpCountsFile));

            // skip field names
            String line = fileReader.readLine();

            if (line == null)
            {
                ISF_LOGGER.error("empty calculated expected rates file({})", mConfig.ExpCountsFile);
                return false;
            }

            int fragLengths = mConfig.FragmentLengthData.size();
            int expectedColCount = ER_COL_CAT + 1 + fragLengths;

            String currentGeneSetId = "";
            String currentTransGeneName = "";
            Map<String, List<CategoryCountsData>> transGeneCategoryData = null;
            List<CategoryCountsData> categoryDataList = null;

            while ((line = fileReader.readLine()) != null)
            {
                String[] items = line.split(",", -1);

                if (items.length != expectedColCount)
                {
                    ISF_LOGGER.error("invalid exp count data length({}) vs expected({}): {}", items.length, expectedColCount, line);
                    mGeneSetCategoryDataMap.clear();
                    return false;
                }

                String geneSetId = items[ER_COL_GENE_SET_ID];
                String transGeneName = items[ER_COL_TRANS_GENE_NAME];
                String categoryStr = items[ER_COL_CAT];

                if(!geneSetId.equals(currentGeneSetId))
                {
                    currentGeneSetId = geneSetId;
                    transGeneCategoryData = Maps.newHashMap();
                    mGeneSetCategoryDataMap.put(geneSetId, transGeneCategoryData);
                    currentTransGeneName = "";
                }

                if(!transGeneName.equals(currentTransGeneName))
                {
                    currentTransGeneName = transGeneName;
                    categoryDataList = Lists.newArrayList();
                    transGeneCategoryData.put(transGeneName, categoryDataList);
                }

                CategoryCountsData catCounts = new CategoryCountsData(categoryStr, fragLengths);
                categoryDataList.add(catCounts);

                for(int i = 0; i < fragLengths; ++i)
                {
                    int count = Integer.parseInt(items[ER_COL_CAT + i + 1]);
                    catCounts.addCounts(count, i);
                }
            }

            ISF_LOGGER.info("loaded {} gene expected counts from file({})",
                    mGeneSetCategoryDataMap.size(), mConfig.ExpCountsFile);
        }
        catch (IOException e)
        {
            ISF_LOGGER.warn("failed to load expected counts file({}): {}", mConfig.ExpCountsFile, e.toString());
            return false;
        }

        return true;
    }


}
