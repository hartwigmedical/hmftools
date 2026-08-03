package com.hartwig.hmftools.compar;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CurationType.EXPECTED;
import static com.hartwig.hmftools.compar.common.CurationType.INVALID;
import static com.hartwig.hmftools.compar.common.CurationType.NONE;
import static com.hartwig.hmftools.compar.common.FileCommon.FLD_CATEGORY;
import static com.hartwig.hmftools.compar.common.FileCommon.FLD_ITEM_KEY;
import static com.hartwig.hmftools.compar.common.FileCommon.FLD_SAMPLE_ID;

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
import com.hartwig.hmftools.compar.common.FieldDisplayInfo;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.KnownMismatch;
import com.hartwig.hmftools.compar.common.MismatchType;

import org.jetbrains.annotations.Nullable;

public class MismatchFile
{
    private static final String FLD_MISMATCH_TYPE = "MismatchType";
    private static final String FLD_OLD_VALUES = "OldValues";
    private static final String FLD_NEW_VALUES = "NewValues";
    private static final String FLD_DIFFERENCES = "Differences";
    private static final String FLD_CURATION_COMMENT = "CurationComment";
    private static final String FLD_CURATION_TYPE = "CurationType";

    // differences is list of the form: field(oldValue/newValue)
    public static String commonHeader(boolean includeSampleId, boolean includeCatagory)
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        if(includeSampleId)
            sj.add(FLD_SAMPLE_ID.toString());

        if(includeCatagory)
            sj.add(FLD_CATEGORY);

        sj.add(FLD_MISMATCH_TYPE).add(FLD_ITEM_KEY).add(FLD_DIFFERENCES);
        return sj.toString();
    }

    public static String header(boolean includeSampleId, boolean includeCuration, boolean includeComment)
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(commonHeader(includeSampleId, true));
        sj.add(FLD_OLD_VALUES);
        sj.add(FLD_NEW_VALUES);

        if(includeCuration)
        {
            sj.add(FLD_CURATION_TYPE);

            if(includeComment)
                sj.add(FLD_CURATION_COMMENT);
        }

        return sj.toString();
    }

    public static String commonTsv(boolean writeCategory, final Mismatch mismatch)
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        if(writeCategory)
        {
            if(mismatch.OldItem != null)
                sj.add(mismatch.OldItem.category().toString());
            else
                sj.add(mismatch.NewItem.category().toString());
        }

        sj.add(mismatch.Type.toString());

        if(mismatch.OldItem != null)
            sj.add(mismatch.OldItem.key());
        else
            sj.add(mismatch.NewItem.key());

        return sj.toString();
    }

    public static String toTsv(final Mismatch mismatch, boolean writeFieldValues, final List<FieldDisplayInfo> fieldDisplayValues)
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        sj.add(commonTsv(!writeFieldValues, mismatch));

        sj.add(String.join(ITEM_DELIM, mismatch.DiffValues));

        if(writeFieldValues)
        {
            for(FieldDisplayInfo fieldInfo : fieldDisplayValues)
            {
                sj.add(fieldInfo.oldValue());
                sj.add(fieldInfo.newValue());
            }
        }
        else
        {
            sj.add(itemValues(mismatch.OldItem, fieldDisplayValues, true));
            sj.add(itemValues(mismatch.NewItem, fieldDisplayValues, false));
        }

        return sj.toString();
    }

    private static String itemValues(final ComparableItem item, final List<FieldDisplayInfo> fieldDisplayValues, boolean useOld)
    {
        if(item == null)
            return "";

        StringJoiner displaySj = new StringJoiner(ITEM_DELIM);

        for(FieldDisplayInfo fieldInfo : fieldDisplayValues)
        {
            displaySj.add(format("%s=%s", fieldInfo.name(), useOld ? fieldInfo.oldValue() : fieldInfo.newValue()));
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

            Integer sampleIndex = fieldsIndexMap.get(FLD_SAMPLE_ID);

            int categoryIndex = fieldsIndexMap.get(FLD_CATEGORY);
            int mismatchTypeIndex = fieldsIndexMap.get(FLD_MISMATCH_TYPE);
            int keyIndex = fieldsIndexMap.get(FLD_ITEM_KEY);
            int differencesIndex = fieldsIndexMap.get(FLD_DIFFERENCES);
            Integer curationIndex = fieldsIndexMap.get(FLD_CURATION_TYPE);
            Integer commentIndex = fieldsIndexMap.get(FLD_CURATION_COMMENT);

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
