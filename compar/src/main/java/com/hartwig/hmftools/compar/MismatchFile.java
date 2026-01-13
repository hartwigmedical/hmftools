package com.hartwig.hmftools.compar;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.DiffFunctions.diffsStr;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.MismatchData;
import com.hartwig.hmftools.compar.common.MismatchType;

import org.jetbrains.annotations.Nullable;

public class MismatchFile
{
    private enum Columns
    {
        SampleId,
        Category,
        MismatchType,
        Key,
        Differences,
        AllValues;
    }

    // Differences is list of the form: field(refValue/otherValue)

    public static String commonHeader(boolean includeSampleId, boolean includeCatagory)
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        if(includeSampleId)
            sj.add(Columns.SampleId.toString());

        if(includeCatagory)
            sj.add(Columns.Category.toString());

        sj.add(Columns.MismatchType.toString()).add(Columns.Key.toString()).add(Columns.Differences.toString());
        return sj.toString();
    }

    public static String header(boolean includeSampleId)
    {
        return commonHeader(includeSampleId, true) + TSV_DELIM + Columns.AllValues;
    }

    public static String commonTsv(boolean writeCategory, final Mismatch mismatch)
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        if(writeCategory)
        {
            if(mismatch.RefItem != null)
                sj.add(mismatch.RefItem.category().toString());
            else
                sj.add(mismatch.NewItem.category().toString());
        }

        sj.add(mismatch.Type.toString());

        if(mismatch.RefItem != null)
            sj.add(mismatch.RefItem.key());
        else
            sj.add(mismatch.NewItem.key());

        return sj.toString();
    }

    public static String toTsv(final Mismatch mismatch, boolean writeFieldValues, final List<String> comparedFieldsNames)
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        sj.add(commonTsv(!writeFieldValues, mismatch));

        sj.add(diffsStr(mismatch.DiffValues));

        if(writeFieldValues)
        {
            final List<String> refFieldValues = mismatch.RefItem != null ? mismatch.RefItem.displayValues() : null;
            final List<String> newFieldValues = mismatch.NewItem != null ? mismatch.NewItem.displayValues() : null;
            int fieldCount = refFieldValues != null ? refFieldValues.size() : newFieldValues.size();

            for(int i = 0; i < fieldCount; ++i)
            {
                if(refFieldValues != null)
                    sj.add(refFieldValues.get(i));
                else
                    sj.add("");

                if(newFieldValues != null)
                    sj.add(newFieldValues.get(i));
                else
                    sj.add("");
            }
        }
        else
        {
            ComparableItem item = mismatch.RefItem != null ? mismatch.RefItem : mismatch.NewItem;

            StringJoiner displaySj = new StringJoiner(ITEM_DELIM);

            List<String> itemDisplayValues = item.displayValues();

            for(int i = 0; i < itemDisplayValues.size(); ++i)
            {
                displaySj.add(format("%s=%s", comparedFieldsNames.get(i), itemDisplayValues.get(i)));
            }

            sj.add(displaySj.toString());
        }

        return sj.toString();
    }

    public static Map<String,List<MismatchData>> loadMismatches(final String filename, @Nullable final String configSampleId)
    {
        Map<String,List<MismatchData>> sampleMismatches = Maps.newHashMap();

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), TSV_DELIM);
            lines.remove(0);

            Integer sampleIndex = fieldsIndexMap.get(Columns.SampleId.toString());

            int categoryIndex = fieldsIndexMap.get(Columns.Category.toString());
            int mismatchTypeIndex = fieldsIndexMap.get(Columns.MismatchType.toString());
            int keyIndex = fieldsIndexMap.get(Columns.Key.toString());
            int differencesIndex = fieldsIndexMap.get(Columns.Differences.toString());

            List<MismatchData> mismatches = null;
            String currentSampleId = "";

            if(sampleIndex == null && configSampleId != null)
            {
                currentSampleId = configSampleId;
                mismatches = Lists.newArrayList();
                sampleMismatches.put(currentSampleId, mismatches);
            }

            for(final String line : lines)
            {
                String[] values = line.split(TSV_DELIM, -1);

                if(sampleIndex != null)
                {
                    String sampleId = values[sampleIndex];

                    if(!sampleId.equals(currentSampleId))
                    {
                        currentSampleId = sampleId;
                        mismatches = sampleMismatches.get(sampleId);

                        if(mismatches == null)
                        {
                            mismatches = Lists.newArrayList();
                            sampleMismatches.put(sampleId, mismatches);
                        }
                    }
                }

                CategoryType categoryType = CategoryType.valueOf( values[categoryIndex]);
                MismatchType mismatchType = MismatchType.valueOf(values[mismatchTypeIndex]);
                String itemKey = values[keyIndex];

                String diffsStr = values[differencesIndex];
                List<String> differences = !diffsStr.isEmpty() ?
                        Arrays.stream(diffsStr.split(ITEM_DELIM, -1)).collect(Collectors.toList()) : Collections.emptyList();

                mismatches.add(new MismatchData(categoryType, mismatchType, itemKey, differences));
            }

            CMP_LOGGER.info("loaded {} mismatches file({})", mismatches.size(), filename);
        }
        catch(IOException e)
        {
            CMP_LOGGER.error("failed to load mismatch file({}): {}", filename, e.toString());
            System.exit(1);
        }

        return sampleMismatches;

    }


}
