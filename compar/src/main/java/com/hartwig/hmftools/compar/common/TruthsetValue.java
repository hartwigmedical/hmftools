package com.hartwig.hmftools.compar.common;

import static java.lang.String.format;

public class TruthsetValue
{
    public final String FieldName;
    public final String Value;
    public final boolean RequireAbsent;

    public TruthsetValue(final String fieldName, final String value, final boolean requireAbsent)
    {
        FieldName = fieldName;
        Value = value;
        RequireAbsent = requireAbsent;
    }

    public String toString() { return format("%s=%s %s", FieldName, Value, RequireAbsent ? "absent" : ""); }
}
