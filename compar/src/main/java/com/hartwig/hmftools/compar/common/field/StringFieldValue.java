package com.hartwig.hmftools.compar.common.field;

import static java.lang.String.format;

public class StringFieldValue extends FieldValue
{
    public final String Value;

    public StringFieldValue(final String name, final String value, final FieldCheck check)
    {
        super(name, check, "%s");
        Value = value;
    }

    @Override
    public boolean checkDifference(final FieldValue other)
    {
        StringFieldValue otherValue = (StringFieldValue)other;

        return !Value.equals(otherValue.Value);
    }

    @Override
    public String toString() { return format("%s=%s", Name, Value); }

    @Override
    public String displayValue() { return format("%s", Value); }
}
