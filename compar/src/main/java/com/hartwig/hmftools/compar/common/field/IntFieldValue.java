package com.hartwig.hmftools.compar.common.field;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

public class IntFieldValue extends FieldValue
{
    public final long Value;

    public IntFieldValue(final FieldInfo field, final long value)
    {
        super(field);
        Value = value;
    }

    @Override
    public boolean checkDifference(final FieldValue other)
    {
        IntFieldValue otherValue = (IntFieldValue)other;

        if(Value == otherValue.Value)
            return false;

        ThresholdFieldCheck thresholdFieldCheck = (ThresholdFieldCheck)Field.fieldCheck();

        double absDiff = abs(Value - otherValue.Value);
        double percDiff = absDiff / max(abs(Value), abs(otherValue.Value));

        boolean satisfiesAbsDiff = thresholdFieldCheck.AbsoluteDiff == null || absDiff > thresholdFieldCheck.AbsoluteDiff;
        boolean satisfiesRelDiff = thresholdFieldCheck.PercentageDiff == null || percDiff > thresholdFieldCheck.PercentageDiff;

        return satisfiesAbsDiff && satisfiesRelDiff;
    }

    @Override
    public String toString() { return format("%s=%d", name(), Value); }

    @Override
    public String displayValue() { return format("%d", Value); }
}
