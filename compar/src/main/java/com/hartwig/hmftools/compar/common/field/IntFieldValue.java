package com.hartwig.hmftools.compar.common.field;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

public class IntFieldValue extends FieldValue
{
    public final long Value;

    public IntFieldValue(final String name, final long value, final ThresholdFieldCheck check)
    {
        super(name, check, "%d");
        Value = value;
    }

    @Override
    public boolean checkDifference(final FieldValue other)
    {
        IntFieldValue otherValue = (IntFieldValue)other;

        if(Value == otherValue.Value)
            return false;

        ThresholdFieldCheck thresholdFieldCheck = (ThresholdFieldCheck)Check;

        double absDiff = abs(Value - otherValue.Value);
        double percDiff = absDiff / max(abs(Value), abs(otherValue.Value));

        boolean satisfiesAbsDiff = thresholdFieldCheck.AbsoluteDiff == null || absDiff > thresholdFieldCheck.AbsoluteDiff;
        boolean satisfiesRelDiff = thresholdFieldCheck.PercentageDiff == null || percDiff > thresholdFieldCheck.PercentageDiff;

        return satisfiesAbsDiff && satisfiesRelDiff;
    }

    @Override
    public String toString() { return format("%s=%d", Name, Value); }

    @Override
    public String displayValue() { return format("%d", Value); }
}
