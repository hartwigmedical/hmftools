package com.hartwig.hmftools.compar.common;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CategoryType.RNA_FUSION;
import static com.hartwig.hmftools.compar.common.CategoryType.SOMATIC_VARIANT;
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

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class ExpectedCache
{
    private final Map<String,Map<CategoryType,List<ExpectedValue>>> mSampleDataMap;

    public static final String FLD_STATUS = "Status";

    private static final String FLD_FIELD_VALUES = "FieldValues";
    private static final String FLD_REMARK = "Remark";

    private static final String FIELD_VALUE_DELIM = ";";
    private static final String FIELD_VALUE_ITEM_DELIM = "=";
    private static final String NOT_ASSESSED_NA = "N/A";

    private static final String VALUE_PRESENT = "Present";
    private static final String VALUE_ABSENT = "Absent";
    private static final String NO_SAMPLE_ID = "";

    // expected format: SampleId (optional),Category,Key,Status (optional),Remark (optional),FieldValues or columns per field value
    // if no row has a SampleId, rows apply to every sample, otherwise only to their own sample

    public ExpectedCache()
    {
        mSampleDataMap = Maps.newHashMap();
    }

    public void loadFiles(final String expectedFilesConfigStr)
    {
        String[] expectedFiles = expectedFilesConfigStr.split(ITEM_DELIM);

        for(String expectedFile : expectedFiles)
        {
            loadExpectedFile(expectedFile);
        }
    }

    public List<ExpectedValue> sampleExpectedEntries(final String sampleId, final CategoryType category)
    {
        return sampleEntries(sampleId, category).stream().filter(x -> !x.RequireAbsent).toList();
    }

    public List<ExpectedValue> sampleAbsentEntries(final String sampleId, final CategoryType category)
    {
        return sampleEntries(sampleId, category).stream().filter(x -> x.RequireAbsent).toList();
    }

    private List<ExpectedValue> sampleEntries(final String sampleId, final CategoryType category)
    {
        Map<CategoryType, List<ExpectedValue>> categoryMap;
        if(mSampleDataMap.size() == 1 && mSampleDataMap.containsKey(NO_SAMPLE_ID))
        {
            categoryMap = mSampleDataMap.get(NO_SAMPLE_ID);
        }
        else
        {
            categoryMap = mSampleDataMap.get(sampleId);
        }
        return categoryMap != null ? categoryMap.getOrDefault(category, Collections.emptyList()) : Collections.emptyList();
    }

    private static final Set<String> FIXED_FIELDS = Set.of(
            FLD_CATEGORY, FLD_SAMPLE_ID, FLD_STATUS, FLD_ITEM_KEY, FLD_FIELD_VALUES, FLD_REMARK);

    private void loadExpectedFile(final String filename)
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

                    Map<CategoryType, List<ExpectedValue>> categoryItemMap = mSampleDataMap.get(sampleId);

                    if(categoryItemMap == null)
                    {
                        categoryItemMap = Maps.newHashMap();
                        mSampleDataMap.put(sampleId, categoryItemMap);
                    }

                    CategoryType categoryType = CategoryType.valueOf(values[categoryIndex]);
                    String itemKey = values[keyIndex];

                    List<ExpectedValue> expectedValues = categoryItemMap.get(categoryType);

                    if(expectedValues == null)
                    {
                        expectedValues = Lists.newArrayList();
                        categoryItemMap.put(categoryType, expectedValues);
                    }

                    boolean requireAbsent = statusIndex != null && values[statusIndex].equalsIgnoreCase(VALUE_ABSENT);

                    if(statusIndex != null)
                        checkStatus(categoryType, values[statusIndex], itemKey);

                    if(fieldValueIndex != null && values[fieldValueIndex].isEmpty())
                    {
                        expectedValues.add(new ExpectedValue(itemKey, FLD_STATUS, "", requireAbsent));
                    }
                    else if(fieldValueIndex != null)
                    {
                        String[] fieldValueItems = values[fieldValueIndex].split(FIELD_VALUE_DELIM, -1);

                        for(String fieldValuePair : fieldValueItems)
                        {
                            String[] svItem = fieldValuePair.split(FIELD_VALUE_ITEM_DELIM, 2);
                            String field = svItem[0];
                            String value = svItem[1];

                            expectedValues.add(new ExpectedValue(itemKey, field, value, requireAbsent));
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
                            expectedValues.add(new ExpectedValue(itemKey, fieldName, value, requireAbsent));
                        }

                        if(!hasFields)
                            expectedValues.add(new ExpectedValue(itemKey, FLD_STATUS, "", requireAbsent));
                    }
                }
                catch(Exception e)
                {
                    CMP_LOGGER.error("invalid expected entry({}) in file({}): {}", line, filename, e.toString());
                    System.exit(1);
                }
            }

            CMP_LOGGER.info("loaded {} expected entries from file({})", lines.size(), filename);
        }
        catch(IOException e)
        {
            CMP_LOGGER.error("failed to load expected file({}): {}", filename, e.toString());
            System.exit(1);
        }
    }

    @VisibleForTesting
    static void checkStatus(final CategoryType categoryType, final String status, final String itemKey)
    {
        boolean requireAbsent = status.equalsIgnoreCase(VALUE_ABSENT);

        if(!requireAbsent && !status.isEmpty() && !status.equalsIgnoreCase(VALUE_PRESENT))
            throw new IllegalArgumentException("Status must be Present, Absent or blank");

        if(requireAbsent && categoryType != SOMATIC_VARIANT && categoryType != RNA_FUSION)
            throw new IllegalArgumentException("Absent is only supported for SOMATIC_VARIANT and RNA_FUSION");

        if(requireAbsent && itemKey.isEmpty())
            throw new IllegalArgumentException("Absent needs a key");

        if(requireAbsent && itemKey.split(":").length > 2)
            throw new IllegalArgumentException("Absent key can't have coordinates");
    }
}
