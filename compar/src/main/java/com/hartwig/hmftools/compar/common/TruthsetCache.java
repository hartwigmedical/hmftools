package com.hartwig.hmftools.compar.common;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.ComparConstants.FLD_CATEGORY;
import static com.hartwig.hmftools.compar.common.ComparConstants.FLD_ITEM_KEY;
import static com.hartwig.hmftools.compar.common.ComparConstants.FLD_SAMPLE_ID;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class TruthsetCache
{
    private final Map<String,Map<CategoryType,List<TruthsetValue>>> mSampleDataMap;

    public static final String FLD_STATUS = "Status";

    private static final String FLD_FIELD_VALUES = "FieldValues";

    private static final String FIELD_VALUE_DELIM = ";";
    private static final String FIELD_VALUE_ITEM_DELIM = "=";
    private static final String NOT_ASSESSED_NA = "N/A";

    private static final String VALUE_PRESENT = "Present";
    private static final String VALUE_ABSENT = "Absent";
    private static final String NO_SAMPLE_ID = "";

    // expected format: SampleId (if has multiple),Category,Key,Status (optional),FieldValues or columns per field value

    public TruthsetCache()
    {
        mSampleDataMap = Maps.newHashMap();
    }

    public void loadFiles(final String truthsetFilesConfigStr)
    {
        String[] truthsetFiles = truthsetFilesConfigStr.split(ITEM_DELIM);

        for(String truthsetFile : truthsetFiles)
        {
            loadTruthsetFile(truthsetFile);
        }
    }

    public List<TruthsetValue> sampleTruthsetEntries(final String sampleId, final CategoryType category)
    {
        Map<CategoryType, List<TruthsetValue>> categoryMap;
        if(mSampleDataMap.size() == 1 && mSampleDataMap.containsKey(NO_SAMPLE_ID))
        {
            categoryMap = mSampleDataMap.get(NO_SAMPLE_ID);
        }
        else
        {
            categoryMap = mSampleDataMap.get(sampleId);
        }
        return categoryMap != null ? categoryMap.get(category) : Collections.emptyList();
    }

    private static final Set<String> FIXED_FIELDS = Set.of(FLD_CATEGORY, FLD_SAMPLE_ID, FLD_STATUS, FLD_ITEM_KEY, FLD_FIELD_VALUES);

    private void loadTruthsetFile(final String filename)
    {
        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), TSV_DELIM);
            lines.remove(0);

            Integer sampleIndex = fieldsIndexMap.get(FLD_SAMPLE_ID);

            int categoryIndex = fieldsIndexMap.get(FLD_CATEGORY);
            int keyIndex = fieldsIndexMap.get(FLD_ITEM_KEY);

            Integer statusIndex = fieldsIndexMap.get(FLD_STATUS);
            Integer fieldValueIndex = fieldsIndexMap.get(FLD_FIELD_VALUES);

            for(final String line : lines)
            {
                String[] values = line.split(TSV_DELIM, -1);

                try
                {
                    String sampleId = sampleIndex != null ? values[sampleIndex] : NO_SAMPLE_ID;

                    Map<CategoryType, List<TruthsetValue>> categoryItemMap = mSampleDataMap.get(sampleId);

                    if(categoryItemMap == null)
                    {
                        categoryItemMap = Maps.newHashMap();
                        mSampleDataMap.put(sampleId, categoryItemMap);
                    }

                    CategoryType categoryType = CategoryType.valueOf(values[categoryIndex]);
                    String itemKey = values[keyIndex];

                    List<TruthsetValue> truthsetValues = categoryItemMap.get(categoryType);

                    if(truthsetValues == null)
                    {
                        truthsetValues = Lists.newArrayList();
                        categoryItemMap.put(categoryType, truthsetValues);
                    }

                    boolean requireAbsent = statusIndex != null && values[statusIndex].equalsIgnoreCase(VALUE_ABSENT);

                    if(fieldValueIndex != null)
                    {
                        String[] fieldValueItems = values[fieldValueIndex].split(FIELD_VALUE_DELIM, -1);

                        for(String fieldValuePair : fieldValueItems)
                        {
                            String[] svItem = fieldValuePair.split(FIELD_VALUE_ITEM_DELIM, 2);
                            String field = svItem[0];
                            String value = svItem[1];

                            truthsetValues.add(new TruthsetValue(itemKey, field, value, requireAbsent));
                        }
                    }
                    else
                    {
                        boolean hasFields = false;

                        for(Map.Entry<String, Integer> entry : fieldsIndexMap.entrySet())
                        {
                            String fieldName = entry.getKey();

                            if(FIXED_FIELDS.contains(fieldName))
                                continue;

                            int fieldIndex = entry.getValue();
                            String value = values[fieldIndex];

                            if(value.equals(NOT_ASSESSED_NA)) // rather than assuming an empty column means not known
                                continue;

                            hasFields = true;
                            truthsetValues.add(new TruthsetValue(itemKey, fieldName, value, requireAbsent));
                        }

                        if(!hasFields)
                            truthsetValues.add(new TruthsetValue(itemKey, FLD_STATUS, "", requireAbsent));
                    }
                }
                catch(Exception e)
                {
                    CMP_LOGGER.error("invalid truthset entry({}}) in file({})", line, filename);
                    System.exit(1);
                }
            }

            CMP_LOGGER.info("loaded {} truthset entries from file({})", lines.size(), filename);
        }
        catch(IOException e)
        {
            CMP_LOGGER.error("failed to load truthset file({}): {}", filename, e.toString());
            System.exit(1);
        }
    }
}
