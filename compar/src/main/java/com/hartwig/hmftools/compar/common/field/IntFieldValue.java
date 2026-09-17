package com.hartwig.hmftools.compar.common.field;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

public class IntFieldValue extends FieldValue
{
    public final int Value;

    public IntFieldValue(final FieldInfo field, final int value)
    {
        super(field);
        Value = value;
    }

    @Override
    public boolean hasDifference(final FieldValue other)
    {
        IntFieldValue otherValue = (IntFieldValue)other;

        if(Value == otherValue.Value)
            return false;

        return Field.FieldCheck.hasThresholds() ? Field.FieldCheck.hasDifference(Value, otherValue.Value) : true;
    }

    @Override
    public String toString() { return format("%s=%d", name(), Value); }

    @Override
    public String displayValue() { return format("%d", Value); }
}
