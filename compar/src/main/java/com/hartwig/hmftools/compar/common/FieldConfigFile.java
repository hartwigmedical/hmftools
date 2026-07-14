package com.hartwig.hmftools.compar.common;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.compar.common.field.DisplayOnlyField.DISPLAY_TYPE;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.compar.common.field.Field;

public class FieldConfigFile
{
    static final String NONE_SETTING = "none";

    private static final String COL_CATEGORY = "Category";
    private static final String COL_FIELD = "Field";
    private static final String COL_FIELD_TYPE = "FieldType";
    public static final String COL_COMPARED = "Compared";
    public static final String COL_ABSOLUTE_THRESHOLD = "AbsoluteThreshold";
    public static final String COL_PERCENT_THRESHOLD = "PercentThreshold";

    public static String generateFileName(final String basePath)
    {
        return checkAddDirSeparator(basePath) + "field.config.compar.tsv";
    }

    public static void write(final String filename, final FieldConfig fieldConfig, final Set<CategoryType> categories) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(fieldConfig, categories));
    }

    public static List<FieldOverride> read(final String filename) throws IOException
    {
        List<String> lines = Files.readAllLines(new File(filename).toPath());

        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), TSV_DELIM);
        lines.remove(0);

        int categoryIndex = fieldsIndexMap.get(COL_CATEGORY);
        int fieldIndex = fieldsIndexMap.get(COL_FIELD);
        int comparedIndex = fieldsIndexMap.get(COL_COMPARED);
        int absoluteThresholdIndex = fieldsIndexMap.get(COL_ABSOLUTE_THRESHOLD);
        int percentThresholdIndex = fieldsIndexMap.get(COL_PERCENT_THRESHOLD);

        List<FieldOverride> fieldOverrides = Lists.newArrayList();

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            fieldOverrides.add(new FieldOverride(
                    values[categoryIndex], values[fieldIndex], values[comparedIndex],
                    values[absoluteThresholdIndex], values[percentThresholdIndex]));
        }

        return fieldOverrides;
    }

    static List<String> toLines(final FieldConfig fieldConfig, final Set<CategoryType> categories)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        for(CategoryType category : CategoryType.values())
        {
            if(categories.contains(category))
            {
                for(Field field: fieldConfig.getFields(category))
                {
                    if(!field.type().equals(DISPLAY_TYPE))
                    {
                        lines.add(toLine(field, category));
                    }
                }
            }
        }
        return lines;
    }

    private static String toLine(final Field field, final CategoryType category)
    {
        return new StringJoiner(TSV_DELIM)
                .add(category.toString())
                .add(field.name())
                .add(field.type())
                .add(String.valueOf(field.isCompared()))
                .add(field.absoluteThreshold() == null ? NONE_SETTING : String.valueOf(field.absoluteThreshold()))
                .add(field.percentThreshold() == null ? NONE_SETTING : (field.percentThreshold() * 100) + "%")
                .toString();
    }

    private static String header()
    {
        return new StringJoiner(TSV_DELIM)
                .add(COL_CATEGORY)
                .add(COL_FIELD)
                .add(COL_FIELD_TYPE)
                .add(COL_COMPARED)
                .add(COL_ABSOLUTE_THRESHOLD)
                .add(COL_PERCENT_THRESHOLD)
                .toString();
    }
}
