package com.hartwig.hmftools.compar.common.field;

import java.util.List;

public class FieldInfo
{
    private final String mName;
    private final String mFormatString;
    private final FieldCheck mFieldCheck;
    private final boolean mDisplayOnly;

    public FieldInfo(final String name, final FieldCheck fieldCheck, final String formatString)
    {
        mName = name;
        mFormatString = formatString;

        if(fieldCheck == null)
        {
            mFieldCheck = new FieldCheck(false);
            mDisplayOnly = true;
        }
        else
        {
            mFieldCheck = fieldCheck;
            mDisplayOnly = false;
        }
    }

    public static FieldInfo displayOnly(final String name, final String formatString)
    {
        return new FieldInfo(name, null, formatString);
    }

    public String name() { return mName; }
    public String formatString() { return mFormatString; }
    public FieldCheck fieldCheck() { return mFieldCheck; }
    public boolean displayOnly() { return mDisplayOnly; }

    public static FieldInfo findField(final String fieldName, final List<FieldInfo> fields)
    {
        return fields.stream().filter(x -> x.name().equals(fieldName)).findFirst().orElse(null);
    }
}
