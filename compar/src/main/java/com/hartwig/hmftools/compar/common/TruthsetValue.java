package com.hartwig.hmftools.compar.common;

import static java.lang.String.format;

public class TruthsetValue
{
    public final String Key;
    public final String FieldName;
    public final String Value;
    public final boolean RequireAbsent;

    public TruthsetValue(final String key, final String fieldName, final String value, final boolean requireAbsent)
    {
        Key = key;
        FieldName = fieldName;
        Value = value;
        RequireAbsent = requireAbsent;
    }

    public String toString() { return format("%s: %s=%s %s", Key, FieldName, Value, RequireAbsent ? "absent" : ""); }
}
