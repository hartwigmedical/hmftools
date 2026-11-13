package com.hartwig.hmftools.compar.common.field;

import static java.lang.String.format;

import java.util.List;

public class FieldInfo
{
    public final String Name;
    public final String FormatString;
    public final FieldCheck FieldCheck;
    public final boolean DisplayOnly;

    public FieldInfo(final String name, final FieldCheck fieldCheck, final String formatString)
    {
        Name = name;
        FormatString = formatString;

        if(fieldCheck == null)
        {
            FieldCheck = new FieldCheck(false);
            DisplayOnly = true;
        }
        else
        {
            FieldCheck = fieldCheck;
            DisplayOnly = false;
        }
    }

    public static FieldInfo displayOnly(final String name, final String formatString)
    {
        return new FieldInfo(name, null, formatString);
    }

    public static FieldInfo findField(final String fieldName, final List<FieldInfo> fields)
    {
        return fields.stream().filter(x -> x.Name.equals(fieldName)).findFirst().orElse(null);
    }

    public String toString() { return format("%s: %s %s", Name, FieldCheck, DisplayOnly ? "display-only" : ""); }
}
