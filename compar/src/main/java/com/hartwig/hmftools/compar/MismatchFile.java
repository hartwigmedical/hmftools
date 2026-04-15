package com.hartwig.hmftools.compar;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CurationType.EXPECTED;
import static com.hartwig.hmftools.compar.common.CurationType.INVALID;
import static com.hartwig.hmftools.compar.common.CurationType.NONE;
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
import com.hartwig.hmftools.compar.common.CurationInfo;
import com.hartwig.hmftools.compar.common.CurationType;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.KnownMismatch;
import com.hartwig.hmftools.compar.common.MismatchType;

import org.jetbrains.annotations.Nullable;

public class MismatchFile
{
    private enum Columns // for the generic file output
    {
        SampleId,
        Category,
        MismatchType,
        Key,
        Differences,
        RefValues,
        NewValues,
        CurationType,
        CurationComment;
    }

    // differences is list of the form: field(refValue/newValue)
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

    public static String header(boolean includeSampleId, boolean includeCuration, boolean includeComment)
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(commonHeader(includeSampleId, true));
        sj.add(Columns.RefValues.toString());
        sj.add(Columns.NewValues.toString());

        if(includeCuration)
        {
            sj.add(Columns.CurationType.toString());

            if(includeComment)
                sj.add(Columns.CurationComment.toString());
        }

        return sj.toString();
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
            List<String> refFieldValues = mismatch.RefItem != null ? mismatch.RefItem.displayValues() : null;
            List<String> newFieldValues = mismatch.NewItem != null ? mismatch.NewItem.displayValues() : null;
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
            sj.add(itemValues(mismatch.RefItem, comparedFieldsNames));
            sj.add(itemValues(mismatch.NewItem, comparedFieldsNames));
        }

        return sj.toString();
    }

    private static String itemValues(final ComparableItem item, final List<String> comparedFieldsNames)
    {
        if(item == null)
            return "";

        StringJoiner displaySj = new StringJoiner(ITEM_DELIM);

        List<String> itemDisplayValues = item.displayValues();

        for(int i = 0; i < itemDisplayValues.size(); ++i)
        {
            displaySj.add(format("%s=%s", comparedFieldsNames.get(i), itemDisplayValues.get(i)));
        }

        return displaySj.toString();
    }

    public static Map<String,List<KnownMismatch>> loadSampleCurations(final String filename, @Nullable final String configSampleId)
    {
        Map<String,List<KnownMismatch>> sampleMismatches = Maps.newHashMap();

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
            Integer curationIndex = fieldsIndexMap.get(Columns.CurationType.toString());
            Integer commentIndex = fieldsIndexMap.get(Columns.CurationComment.toString());

            List<KnownMismatch> mismatches = null;
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

                CategoryType categoryType = CategoryType.valueOf(values[categoryIndex]);
                MismatchType mismatchType = MismatchType.valueOf(values[mismatchTypeIndex]);
                String itemKey = values[keyIndex];

                String diffsStr = values[differencesIndex];

                CurationType curationType = NONE;
                String curationComment = "";

                if(curationIndex != null)
                {
                    try
                    {
                        if(curationIndex < values.length)
                            curationType = CurationType.valueOf(values[curationIndex]);

                        if(commentIndex != null && commentIndex < values.length)
                            curationComment = values[commentIndex];
                    }
                    catch(Exception e)
                    {
                        curationType = INVALID; // use has configured an invalid enum/type
                    }
                }
                else
                {
                    curationType = EXPECTED; // loading a prior unedited mismatch file
                }

                List<String> differences = !diffsStr.isEmpty() ?
                        Arrays.stream(diffsStr.split(ITEM_DELIM, -1)).collect(Collectors.toList()) : Collections.emptyList();

                CurationInfo curationInfo = new CurationInfo(curationType, curationComment);

                mismatches.add(new KnownMismatch(categoryType, mismatchType, itemKey, differences, curationInfo));
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
