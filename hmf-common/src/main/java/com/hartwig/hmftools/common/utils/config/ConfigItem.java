package com.hartwig.hmftools.common.utils.config;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.config.ConfigItemType.FLAG;

public class ConfigItem
{
    public final ConfigItemType Type;
    public final String Name;
    public final boolean Required;
    public final String Description;

    private String mValue;
    private boolean mHasValue;
    private String mDefaultValue;

    public ConfigItem(final ConfigItemType type, final String name, final boolean required, final String description,
            final String defaultValue)
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

        mValue = mDefaultValue != null ? mDefaultValue : "";
        mHasValue = false;
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

    public boolean missing() { return Required && mValue.isEmpty(); }

    public void clearValue()
    {
        mValue = mDefaultValue;
        mHasValue = false;
    }

    public String toString()
    {
        return format("%s:%s required(%s) value(%s) default(%s) desc(%s)",
                Type, Name, Required, mValue, mDefaultValue != null ? mDefaultValue : "none", Description);
    }
}
