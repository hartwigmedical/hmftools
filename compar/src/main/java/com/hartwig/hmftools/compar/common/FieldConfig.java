package com.hartwig.hmftools.compar.common;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.FieldConfigFile.NONE_SETTING;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.field.Field;

public class FieldConfig
{
    private final Map<CategoryType, Map<String, Field>> fieldSettings;

    public FieldConfig()
    {
        fieldSettings = Maps.newHashMap();
    }

    public void registerFields(final ItemComparer comparer, final MatchLevel matchLevel)
    {
        for(Field field : comparer.fields(matchLevel))
        {
            registerField(comparer.category(), field);
        }
    }

    public void registerField(final CategoryType category, final Field field)
    {
        fieldSettings.putIfAbsent(category, Maps.newHashMap());
        fieldSettings.get(category).put(field.name(), field);
    }

    public List<Field> getFields(final CategoryType category)
    {
        return fieldSettings.getOrDefault(category, Maps.newHashMap()).values().stream().toList();
    }

    public List<Field> getFields(final CategoryType category, final List<String> fieldNames)
    {
        return fieldNames.stream().map(n -> fieldSettings.get(category).get(n)).toList();
    }

    public void applyOverride(final FieldOverride fieldOverride)
    {
        CategoryType category;

        try
        {
            category = CategoryType.valueOf(fieldOverride.Category);
        }
        catch(IllegalArgumentException e)
        {
            CMP_LOGGER.warn("field override category({}) does not exist", fieldOverride.Category);
            return;
        }

        Field field = fieldSettings.getOrDefault(category, Maps.newHashMap()).get(fieldOverride.Field);

        if(field == null)
        {
            CMP_LOGGER.warn("field override category({}) field({}) not registered", fieldOverride.Category, fieldOverride.Field);
            return;
        }

        if(!fieldOverride.Compared.isEmpty())
        {
            field = field.withCompared(Boolean.parseBoolean(fieldOverride.Compared));
        }

        if(!fieldOverride.AbsoluteThreshold.isEmpty())
        {
            field = field.withAbsoluteThreshold(parseThreshold(fieldOverride.AbsoluteThreshold));
        }

        if(!fieldOverride.PercentThreshold.isEmpty())
        {
            field = field.withPercentThreshold(parsePercentThreshold(fieldOverride.PercentThreshold));
        }

        fieldSettings.get(category).put(field.name(), field);
    }

    private static Double parseThreshold(final String value)
    {
        return value.equals(NONE_SETTING) ? null : Double.parseDouble(value);
    }

    private static Double parsePercentThreshold(final String value)
    {
        if(value.equals(NONE_SETTING))
        {
            return null;
        }

        if(value.endsWith("%"))
        {
            return Double.parseDouble(value.substring(0, value.length() - 1)) / 100;
        }

        return Double.parseDouble(value);
    }
}
