package com.hartwig.hmftools.compar.common.field;

import static java.lang.String.format;

public class BoolFieldValue extends FieldValue
{
    public final double Value;

    public BoolFieldValue(final FieldInfo field, final double value)
    {
        super(field);
        Value = value;
    }

    @Override
    public boolean checkDifference(final FieldValue other)
    {
        BoolFieldValue otherValue = (BoolFieldValue)other;

        return Value != otherValue.Value;
    }

    @Override
    public String toString() { return format("%s=%s", name(), formatString()); }

    @Override
    public String displayValue() { return format("%s", Value); }
}
