package com.hartwig.hmftools.common.utils.config;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.config.ConfigItemType.FLAG;

import java.util.Arrays;
import java.util.stream.Collectors;

public class ConfigItem
{
    public final ConfigItemType Type;
    public final String Name;
    public final boolean Required;
    public final String Description;

    private String mValue;
    private boolean mHasValue;
    private String mDefaultValue;
    private String mPathPrefixConfig;

    public ConfigItem(
            final ConfigItemType type, final String name, final boolean required, final String description, final String defaultValue)
    {
        Type = type;
        Name = name;
        Description = description;

        if(Type == FLAG)
        {
            Required = false;
            mDefaultValue = Boolean.FALSE.toString();
        }
        else
        {
            Required = required;
            mDefaultValue = defaultValue;
        }

        mValue = mDefaultValue; // can be null
        mHasValue = false;
        mPathPrefixConfig = null;
    }

    public void setValue(final String value)
    {
        mValue = value;
        mHasValue = true;
    }

    public String value() { return mValue; }
    public boolean hasValue() { return mHasValue; }
    public String defaultValue() { return mDefaultValue; }

    public double decimal() { return Double.parseDouble(mValue); }
    public int integer() { return Integer.parseInt(mValue); }
    public boolean bool() { return Boolean.parseBoolean(mValue); }

    public boolean missing() { return Required && !mHasValue; }

    public void clearValue()
    {
        mValue = mDefaultValue;
        mHasValue = false;
    }

    public String pathPrefixName() { return mPathPrefixConfig; }
    public void setPathPrefixName(final String prefix) { mPathPrefixConfig = prefix; }

    public String toString()
    {
        return format("%s %s required(%s) value(%s) default(%s) desc(%s)",
                Name, Type, Required, mValue, mDefaultValue != null ? mDefaultValue : "none", Description);
    }

    public static String enumValuesAsStr(final Enum<?>[] enumValues, final String configName, final String configDelim)
    {
        String valuesStr = Arrays.stream(enumValues).map(x -> x.toString()).collect(Collectors.joining(configDelim));
        return format("%s: %s separated by '%s'", configName, valuesStr, configDelim);
    }
}
