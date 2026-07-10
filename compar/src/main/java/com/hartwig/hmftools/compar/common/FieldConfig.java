package com.hartwig.hmftools.compar.common;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.FieldConfigFile.COL_ABSOLUTE_THRESHOLD;
import static com.hartwig.hmftools.compar.common.FieldConfigFile.COL_COMPARED;
import static com.hartwig.hmftools.compar.common.FieldConfigFile.COL_PERCENT_THRESHOLD;
import static com.hartwig.hmftools.compar.common.FieldConfigFile.NONE_SETTING;
import static com.hartwig.hmftools.compar.common.field.DisplayField.DISPLAY_TYPE;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.field.Field;
import com.hartwig.hmftools.compar.common.field.UnsupportedFieldOverrideException;

public class FieldConfig
{
    private final Map<CategoryType, Map<String, Field>> fieldSettings;
    private final List<String> mWarnings;
    private final List<String> mErrorMessages;

    public FieldConfig()
    {
        fieldSettings = Maps.newHashMap();
        mWarnings = Lists.newArrayList();
        mErrorMessages = Lists.newArrayList();
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

    public void applyOverrides(final List<FieldOverride> overrides, final boolean strictFieldConfig)
    {
        for(FieldOverride fieldOverride : overrides)
        {
            applyOverride(fieldOverride, strictFieldConfig);
        }

        if(strictFieldConfig)
        {
            checkAllFieldsOverridden(overrides);
        }
    }

    public boolean hasErrors() { return !mErrorMessages.isEmpty(); }

    public void logProblems()
    {
        mWarnings.forEach(CMP_LOGGER::warn);
        mErrorMessages.forEach(CMP_LOGGER::error);
    }

    List<String> warnings() { return mWarnings; }
    List<String> errorMessages() { return mErrorMessages; }

    private void applyOverride(final FieldOverride fieldOverride, final boolean strictFieldConfig)
    {
        CategoryType category;

        try
        {
            category = CategoryType.valueOf(fieldOverride.Category);
        }
        catch(IllegalArgumentException e)
        {
            recordProblem(String.format("field override for unknown: category(%s)", fieldOverride.Category), strictFieldConfig);
            return;
        }

        Field field = fieldSettings.getOrDefault(category, Maps.newHashMap()).get(fieldOverride.Field);

        if(field == null)
        {
            String message =
                    String.format("field override for unknown field: category(%s) field(%s)", fieldOverride.Category, fieldOverride.Field);
            recordProblem(message, strictFieldConfig);
            return;
        }

        field = applyComparedOverride(fieldOverride, field, category, strictFieldConfig);
        field = applyAbsoluteThresholdOverride(fieldOverride, field, category, strictFieldConfig);
        field = applyPercentThresholdOverride(fieldOverride, field, category, strictFieldConfig);

        fieldSettings.get(category).put(field.name(), field);
    }

    private Field applyComparedOverride(
            final FieldOverride fieldOverride, final Field field, final CategoryType category, final boolean strictFieldConfig)
    {
        if(!fieldOverride.Compared.isEmpty())
        {
            Boolean compared = parseCompared(fieldOverride.Compared, field.name());

            if(compared != null)
            {
                try
                {
                    return field.withCompared(compared);
                }
                catch(UnsupportedFieldOverrideException e)
                {
                    mErrorMessages.add(e.getMessage());
                }
            }
        }
        else if(strictFieldConfig)
        {
            recordMissingOverrideError(field, category, COL_COMPARED);
        }

        return field;
    }

    private Field applyAbsoluteThresholdOverride(
            final FieldOverride fieldOverride, final Field field, final CategoryType category, final boolean strictFieldConfig)
    {
        if(!fieldOverride.AbsoluteThreshold.isEmpty())
        {
            Double absoluteThreshold = parseAbsoluteThreshold(fieldOverride.AbsoluteThreshold, field.name());

            try
            {
                return field.withAbsoluteThreshold(absoluteThreshold);
            }
            catch(UnsupportedFieldOverrideException e)
            {
                mErrorMessages.add(e.getMessage());
            }
        }
        else if(strictFieldConfig)
        {
            recordMissingOverrideError(field, category, COL_ABSOLUTE_THRESHOLD);
        }

        return field;
    }

    private Field applyPercentThresholdOverride(
            final FieldOverride fieldOverride, final Field field, final CategoryType category, final boolean strictFieldConfig)
    {
        if(!fieldOverride.PercentThreshold.isEmpty())
        {
            Double percentThreshold = parsePercentThreshold(fieldOverride.PercentThreshold, field.name());

            try
            {
                return field.withPercentThreshold(percentThreshold);
            }
            catch(UnsupportedFieldOverrideException e)
            {
                mErrorMessages.add(e.getMessage());
            }
        }
        else if(strictFieldConfig)
        {
            recordMissingOverrideError(field, category, COL_PERCENT_THRESHOLD);
        }

        return field;
    }

    private void recordMissingOverrideError(final Field field, final CategoryType category, final String columnName)
    {
        mErrorMessages.add(String.format(
                "field(%s) category(%s) missing '%s' override in strict field config mode", field.name(), category, columnName));
    }

    private void checkAllFieldsOverridden(final List<FieldOverride> overrides)
    {
        Set<String> overriddenKeys = overrides.stream()
                .map(o -> o.Category + ":" + o.Field)
                .collect(Collectors.toSet());

        for(Map.Entry<CategoryType, Map<String, Field>> entry : fieldSettings.entrySet())
        {
            CategoryType category = entry.getKey();

            for(Field field : entry.getValue().values())
            {
                if(field.type().equals(DISPLAY_TYPE))
                {
                    continue;
                }

                String key = category + ":" + field.name();

                if(!overriddenKeys.contains(key))
                {
                    mErrorMessages.add(String.format(
                            "field config file missing entry for category(%s) field(%s)", category, field.name()));
                }
            }
        }
    }

    private void recordProblem(final String message, final boolean isError)
    {
        if(isError)
        {
            mErrorMessages.add(message);
        }
        else
        {
            mWarnings.add(message);
        }
    }

    private Boolean parseCompared(final String value, final String fieldName)
    {
        if(value.equalsIgnoreCase("true"))
        {
            return Boolean.TRUE;
        }

        if(value.equalsIgnoreCase("false"))
        {
            return Boolean.FALSE;
        }

        mErrorMessages.add(String.format("compared override for field %s cannot be parsed: %s", fieldName, value));
        return null;
    }

    private Double parseAbsoluteThreshold(final String value, final String fieldName)
    {
        Double result;
        try{
            if(value.equals(NONE_SETTING))
            {
                result = null;
            }
            else
            {
                result = Double.parseDouble(value);
            }
        }
        catch(NumberFormatException e)
        {
            mErrorMessages.add(String.format("absolute threshold for field %s cannot be parsed: %s", fieldName, value));
            result = null;
        }

        return result;
    }

    private Double parsePercentThreshold(final String value, final String fieldName)
    {
        Double result;
        try{
            if(value.equals(NONE_SETTING))
            {
                result = null;
            }
            else if(value.endsWith("%"))
            {
                result = Double.parseDouble(value.substring(0, value.length() - 1)) / 100;
            }
            else
            {
                result = Double.parseDouble(value);
            }
        }
        catch(NumberFormatException e)
        {
            mErrorMessages.add(String.format("percent threshold for field %s cannot be parsed: %s", fieldName, value));
            result = null;
        }

        return result;
    }
}
