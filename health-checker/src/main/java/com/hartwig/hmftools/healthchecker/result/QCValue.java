package com.hartwig.hmftools.healthchecker.result;

public class QCValue
{
    public final QCValueType Type;
    public final String Value;

    public QCValue(final QCValueType type, final String value)
    {
        Type = type;
        Value = value;
    }
}
