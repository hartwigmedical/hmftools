package com.hartwig.hmftools.compar.common.field;

import static java.lang.String.format;

public class StringFieldValue extends FieldValue
{
    public final String Value;

    public StringFieldValue(final FieldInfo field, final String value)
    {
        super(field);
        Value = value;
    }

    @Override
    public boolean hasDifference(final FieldValue other)
    {
        StringFieldValue otherValue = (StringFieldValue)other;

        return !Value.equals(otherValue.Value);
    }

    @Override
    public String toString() { return format("%s=%s", name(), Value); }

    @Override
    public String displayValue() { return Value; }
}
