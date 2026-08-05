package com.hartwig.hmftools.compar.common.field;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

public class DoubleFieldValue extends FieldValue
{
    public final double Value;

    public DoubleFieldValue(final FieldInfo field, final double value)
    {
        super(field);
        Value = value;
    }

    @Override
    public boolean hasDifference(final FieldValue other)
    {
        DoubleFieldValue otherDouble = (DoubleFieldValue)other;

        if(Value == otherDouble.Value)
            return false;

        ThresholdFieldCheck thresholdFieldCheck = (ThresholdFieldCheck)Field.FieldCheck;

        double absDiff = abs(Value - otherDouble.Value);
        double percDiff = absDiff / max(abs(Value), abs(otherDouble.Value));

        boolean satisfiesAbsDiff = thresholdFieldCheck.AbsoluteDiff == null || absDiff > thresholdFieldCheck.AbsoluteDiff;
        boolean satisfiesRelDiff = thresholdFieldCheck.PercentageDiff == null || percDiff > thresholdFieldCheck.PercentageDiff;

        return satisfiesAbsDiff && satisfiesRelDiff;
    }

    @Override
    public String toString() { return format("%s=%s", name(), format(formatString(), Value)); }

    @Override
    public String displayValue() { return format("%s", format(formatString(), Value)); }
}
