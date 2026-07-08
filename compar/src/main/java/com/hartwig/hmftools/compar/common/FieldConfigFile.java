package com.hartwig.hmftools.compar.common;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.compar.common.field.DisplayField.DISPLAY_TYPE;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.compar.common.field.Field;

public class FieldConfigFile
{
    private static final String NONE_SETTING = "none";

    public static String generateFileName(final String basePath)
    {
        return checkAddDirSeparator(basePath) + "field.config.compar.tsv";
    }

    public static void write(final String filename, final FieldConfig fieldConfig, final Set<CategoryType> categories) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(fieldConfig, categories));
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
                .add("category")
                .add("field")
                .add("fieldType")
                .add("compared")
                .add("absoluteThreshold")
                .add("percentThreshold")
                .toString();
    }
}
