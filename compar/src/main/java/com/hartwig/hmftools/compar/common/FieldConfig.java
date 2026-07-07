package com.hartwig.hmftools.compar.common;

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
}
