package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.FileCommon.FLD_CATEGORY;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.field.FieldInfo;
import com.hartwig.hmftools.compar.common.field.ThresholdFieldCheck;

public class FieldConfigFile
{
    static final String NONE_SETTING = "none";

    public static final String FLD_FIELD = "Field";
    public static final String FLD_COMPARED = "Compared";
    public static final String FLD_ABSOLUTE_THRESHOLD = "AbsoluteThreshold";
    public static final String FLD_PERCENT_THRESHOLD = "PercentThreshold";

    public static String generateFileName(final String filePrefix)
    {
        return filePrefix + ".field_config.tsv";
    }

    public static void write(final String filename, List<ItemComparer> comparers)
    {
        try
        {
            List<String> lines = Lists.newArrayList();
            lines.add(header());

            for(ItemComparer comparer : comparers)
            {
                List<FieldInfo> fields = comparer.fieldsList();

                fields.stream().filter(x -> !x.DisplayOnly).forEach(x -> lines.add(toLine(comparer.category(), x)));
            }

            Files.write(new File(filename).toPath(), lines);
        }
        catch(IOException e)
        {
            CMP_LOGGER.error("error writing field config file", e.toString());
            System.exit(1);
        }
    }

    private static String toLine(final CategoryType category, final FieldInfo fieldInfo)
    {
        Double absoluteThreshold = null;
        Double percentThreshold = null;

        if(fieldInfo.FieldCheck instanceof ThresholdFieldCheck)
        {
            ThresholdFieldCheck thresholdFieldCheck = (ThresholdFieldCheck)fieldInfo.FieldCheck;
            absoluteThreshold = thresholdFieldCheck.AbsoluteDiff;
            percentThreshold = thresholdFieldCheck.PercentageDiff;
        }

        return new StringJoiner(TSV_DELIM)
                .add(category.toString())
                .add(fieldInfo.Name)
                .add(String.valueOf(fieldInfo.FieldCheck.IsCompared))
                .add(absoluteThreshold == null ? NONE_SETTING : String.valueOf(absoluteThreshold))
                .add(percentThreshold == null ? NONE_SETTING : (percentThreshold * 100) + "%")
                .toString();
    }

    private static String header()
    {
        return new StringJoiner(TSV_DELIM)
                .add(FLD_CATEGORY)
                .add(FLD_FIELD)
                .add(FLD_COMPARED)
                .add(FLD_ABSOLUTE_THRESHOLD)
                .add(FLD_PERCENT_THRESHOLD)
                .toString();
    }
}
