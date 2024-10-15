package com.hartwig.hmftools.healthchecker;

import static java.lang.String.format;

public class QCValue
{
    public final QCValueType Type;
    public final String Value;

    public QCValue(final QCValueType type, final String value)
    {
        Type = type;
        Value = value;
    }

    public String toString() { return format("%s=%s", Type, Value); }
}
