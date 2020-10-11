package com.hartwig.hmftools.isofox.expression;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.expression.ExpectedRatesGenerator.EXP_COUNT_LENGTH_HEADER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

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
import com.hartwig.hmftools.isofox.adjusts.FragmentSize;

public class ExpectedCountsCache
{
    private final IsofoxConfig mConfig;
    private final Map<String,Map<String,List<CategoryCountsData>>> mGeneSetCategoryDataMap;
    private boolean mValidData;

    public ExpectedCountsCache(final IsofoxConfig config)
    {
        mConfig = config;
        mGeneSetCategoryDataMap = Maps.newHashMap();
        mValidData = true;

        if(config.ExpCountsFile != null && Files.exists(Paths.get(mConfig.ExpCountsFile)))
        {
            mValidData = loadExpCountsFile();
        }
    }

    public Map<String,List<CategoryCountsData>> getGeneExpectedRatesData(final String chrId, final List<String> geneIds)
    {
        Map<String,List<CategoryCountsData>> geneSetCountsData = mGeneSetCategoryDataMap.get(chrId);

        if (geneSetCountsData == null || !geneSetCountsDataMatches(geneIds, geneSetCountsData.keySet()))
        {
            geneSetCountsData = findGeneSetCountsData(geneIds);
        }

        return geneSetCountsData;
    }

    private boolean geneSetCountsDataMatches(final List<String> geneIds, final Set<String> geneTransSet)
    {
        // confirm that the genes in the collection match
        return !geneIds.stream().anyMatch(x -> !geneTransSet.contains(x));
    }

    private final Map<String,List<CategoryCountsData>> findGeneSetCountsData(final List<String> geneIds)
    {
        for(final Map<String,List<CategoryCountsData>> geneCounts : mGeneSetCategoryDataMap.values())
        {
            if(geneSetCountsDataMatches(geneIds, geneCounts.keySet()))
                return geneCounts;
        }

        return null;
    }

    // GeneSetId,TransId,Category,Counts for each fragment length
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
                ISF_LOGGER.error("empty calculated expected counts file({})", mConfig.ExpCountsFile);
                return false;
            }

            String[] headerItems = line.split(DELIMITER, -1);

            // extract the fragment lengths from the header if not already populated (in which case they must match)
            int fileFragmentLengthCount = headerItems.length - 3;

            if(mConfig.FragmentSizeData.size() == 0)
            {
                for(int i = 3; i < headerItems.length; ++i)
                {
                    int fragmentLength = Integer.parseInt(headerItems[i].replaceAll(EXP_COUNT_LENGTH_HEADER, ""));
                    mConfig.FragmentSizeData.add(new FragmentSize(fragmentLength, 0));
                }
            }
            else if(mConfig.FragmentSizeData.size() != fileFragmentLengthCount)
            {
                ISF_LOGGER.error("expected counts file has {} fragment lengths vs configuredCount({})",
                        fileFragmentLengthCount, mConfig.FragmentSizeData.size());
                return false;
            }

            int fragLengths = mConfig.FragmentSizeData.size();

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(line, DELIMITER);
            int geneSetIdIndex = fieldsIndexMap.get("GeneSetId");
            int transNameIndex = fieldsIndexMap.get("TransId");
            int categoryIndex = fieldsIndexMap.get("Category");

            String currentGeneSetId = "";
            String currentTransGeneName = "";
            Map<String, List<CategoryCountsData>> transGeneCategoryData = null;
            List<CategoryCountsData> categoryDataList = null;

            while ((line = fileReader.readLine()) != null)
            {
                String[] items = line.split(DELIMITER, -1);

                String geneSetId = items[geneSetIdIndex];
                String transGeneName = items[transNameIndex];
                String categoryStr = items[categoryIndex];

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
                    int count = Integer.parseInt(items[categoryIndex + i + 1]);
                    catCounts.addFragLengthCounts(count, i);
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
