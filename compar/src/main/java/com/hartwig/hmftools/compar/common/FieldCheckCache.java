package com.hartwig.hmftools.compar.common;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.FieldConfigFile.FLD_ABSOLUTE_THRESHOLD;
import static com.hartwig.hmftools.compar.common.FieldConfigFile.FLD_COMPARED;
import static com.hartwig.hmftools.compar.common.FieldConfigFile.FLD_FIELD;
import static com.hartwig.hmftools.compar.common.FieldConfigFile.FLD_FIELD_TYPE;
import static com.hartwig.hmftools.compar.common.FieldConfigFile.FLD_PERCENT_THRESHOLD;
import static com.hartwig.hmftools.compar.common.FieldConfigFile.NONE_SETTING;
import static com.hartwig.hmftools.compar.common.FileCommon.FLD_CATEGORY;
import static com.hartwig.hmftools.compar.common.field.DisplayOnlyField.DISPLAY_TYPE;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.field.Field;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.ThresholdFieldCheck;
import com.hartwig.hmftools.compar.common.field.UnsupportedFieldOverrideException;
import com.hartwig.hmftools.compar.purple.PurityComparer;

public class FieldCheckCache
{
    private final Map<CategoryType, Map<String, Field>> mFieldSettings;

    private final Map<CategoryType,Map<String,FieldCheck>> mFieldCheckOverrides;

    private final List<String> mWarnings;
    private final List<String> mErrorMessages;

    public FieldCheckCache()
    {
        // use linked hash map to preserve field registration order in diffs
        mFieldSettings = Maps.newLinkedHashMap();
        mFieldCheckOverrides = Maps.newHashMap();
        mWarnings = Lists.newArrayList();
        mErrorMessages = Lists.newArrayList();
    }

    public Map<String,FieldCheck> getCategoryFieldCheckOverrides(final CategoryType category)
    {
        return mFieldCheckOverrides.containsKey(category) ? mFieldCheckOverrides.get(category) : Collections.emptyMap();
    }

    public List<String> warnings() { return mWarnings; }
    public List<String> errorMessages() { return mErrorMessages; }

    private static final String NOT_ASSESSED_NA = "N/A";

    public void loadOverrides(final String filename) throws IOException
    {
        List<String> lines = Files.readAllLines(new File(filename).toPath());

        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), TSV_DELIM);
        lines.remove(0);

        int categoryIndex = fieldsIndexMap.get(FLD_CATEGORY);
        int fieldIndex = fieldsIndexMap.get(FLD_FIELD);
        int comparedIndex = fieldsIndexMap.get(FLD_COMPARED);
        int absoluteThresholdIndex = fieldsIndexMap.get(FLD_ABSOLUTE_THRESHOLD);
        int percentThresholdIndex = fieldsIndexMap.get(FLD_PERCENT_THRESHOLD);

        List<FieldOverride> fieldOverrides = Lists.newArrayList();

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            try
            {
                CategoryType categoryType = CategoryType.valueOf(values[categoryIndex]);
                String fieldName = values[fieldIndex];
                boolean isCompared = Boolean.parseBoolean(values[comparedIndex]);

                Double absThreshold = parseThreshold(values[absoluteThresholdIndex]);
                Double percThreshold = parseThreshold(values[percentThresholdIndex]);

                FieldCheck fieldCheck;

                if(absThreshold != null || percThreshold != null)
                {
                    fieldCheck = new ThresholdFieldCheck(isCompared, absThreshold, percThreshold);
                }
                else
                {
                    fieldCheck = new FieldCheck(isCompared);
                }

                Map<String,FieldCheck> categoryOverrides = mFieldCheckOverrides.get(categoryType);

                if(categoryOverrides == null)
                {
                    categoryOverrides = Maps.newHashMap();
                    mFieldCheckOverrides.put(categoryType, categoryOverrides);
                }

                categoryOverrides.put(fieldName, fieldCheck);
            }
            catch(Exception e)
            {
                CMP_LOGGER.error("invalid field override entry: {}", line);
                mErrorMessages.add(e.getMessage());
            }
        }

        // TODO: handle requirement for 'strict' / comprehensive field overrides

        CMP_LOGGER.info("loaded {} field check overrides from {}", lines.size(), filename);
    }

    private static Double parseThreshold(final String value)
    {
        if(value.isEmpty() || value.equals(NOT_ASSESSED_NA))
            return null;

        return Double.parseDouble(value);
    }

    /* TODO: remove
    public void registerFields(final ItemComparer comparer, final MatchLevel matchLevel)
    {
        for(Field field : comparer.fields(matchLevel))
        {
            registerField(comparer.category(), field);
        }
    }
    */

    public void registerField(final CategoryType category, final Field field)
    {
        mFieldSettings.putIfAbsent(category, Maps.newLinkedHashMap());
        mFieldSettings.get(category).put(field.name(), field);
    }

    public List<Field> getFields(final CategoryType category)
    {
        return mFieldSettings.getOrDefault(category, Maps.newLinkedHashMap()).values().stream().toList();
    }

    public List<Field> getFields(final CategoryType category, final List<String> fieldNames)
    {
        return fieldNames.stream().map(n -> mFieldSettings.get(category).get(n)).toList();
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

    public static ThresholdFieldCheck getOrMakeFieldCheck(
            final Map<String,FieldCheck> fieldCheckMap, final String field, final Double absThreshold, final Double percThreshold)
    {
        FieldCheck override = fieldCheckMap.get(field);
        return override != null ? (ThresholdFieldCheck)override : new ThresholdFieldCheck(true, absThreshold, percThreshold);
    }

    public static FieldCheck getOrMakeFieldCheck(final Map<String,FieldCheck> fieldCheckMap, final String field)
    {
        FieldCheck override = fieldCheckMap.get(field);
        return override != null ? override : new FieldCheck(true);
    }

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

        Field field = mFieldSettings.getOrDefault(category, Maps.newHashMap()).get(fieldOverride.Field);

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

        mFieldSettings.get(category).put(field.name(), field);
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
            recordMissingOverrideError(field, category, FLD_COMPARED);
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
            recordMissingOverrideError(field, category, FLD_ABSOLUTE_THRESHOLD);
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
            recordMissingOverrideError(field, category, FLD_PERCENT_THRESHOLD);
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

        for(Map.Entry<CategoryType, Map<String, Field>> entry : mFieldSettings.entrySet())
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
