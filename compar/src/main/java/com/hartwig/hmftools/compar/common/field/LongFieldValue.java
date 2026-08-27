package com.hartwig.hmftools.compar.common.field;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

public class LongFieldValue extends FieldValue
{
    public final long Value;

    public LongFieldValue(final FieldInfo field, final long value)
    {
        super(field);
        Value = value;
    }

    @Override
    public boolean hasDifference(final FieldValue other)
    {
        LongFieldValue otherValue = (LongFieldValue)other;

        if(Value == otherValue.Value)
            return false;

        return Field.FieldCheck.hasThresholds() ? Field.FieldCheck.hasDifference(Value, otherValue.Value) : true;
    }

    @Override
    public String toString() { return format("%s=%d", name(), Value); }

    @Override
    public String displayValue() { return format("%d", Value); }
}
