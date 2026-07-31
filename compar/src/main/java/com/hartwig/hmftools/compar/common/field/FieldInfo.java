package com.hartwig.hmftools.compar.common.field;

import java.util.List;

public class FieldInfo
{
    private final String mName;
    private final String mFormatString;
    private FieldCheck mFieldCheck;

    public FieldInfo(final String name, final FieldCheck check, final String formatString)
    {
        mName = name;
        mFieldCheck = check;
        mFormatString = formatString;
    }

    public String name() { return mName; }
    public String formatString() { return mFormatString; }
    public FieldCheck fieldCheck() { return mFieldCheck; }

    public void overrideCheck(final FieldCheck fieldCheck) { mFieldCheck = fieldCheck; }

    public static FieldInfo findField(final String fieldName, final List<FieldInfo> fields)
    {
        return fields.stream().filter(x -> x.name().equals(fieldName)).findFirst().orElse(null);
    }
}
