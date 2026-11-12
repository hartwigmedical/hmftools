package com.hartwig.hmftools.compar.common;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;

import java.util.List;

import com.google.common.collect.Lists;

public class ExpectedValue
{
    public final String Key;
    public final String FieldName;
    public final String Value;
    public final boolean RequireAbsent;
    public final String Remark; // the expected row's notes, never compared

    public ExpectedValue(final String key, final String fieldName, final String value, final boolean requireAbsent)
    {
        this(key, fieldName, value, requireAbsent, "");
    }

    public ExpectedValue(
            final String key, final String fieldName, final String value, final boolean requireAbsent, final String remark)
    {
        Key = key;
        FieldName = fieldName;
        Value = value;
        RequireAbsent = requireAbsent;
        Remark = remark;
    }

    // the remarks of an item's expected rows, in file order and without repeats
    public static String remark(final List<ExpectedValue> expectedValues)
    {
        List<String> remarks = Lists.newArrayList();

        for(ExpectedValue expectedValue : expectedValues)
        {
            if(!expectedValue.Remark.isEmpty() && !remarks.contains(expectedValue.Remark))
                remarks.add(expectedValue.Remark);
        }

        return String.join(ITEM_DELIM, remarks);
    }

    public String toString() { return format("%s: %s=%s %s", Key, FieldName, Value, RequireAbsent ? "absent" : ""); }
}
