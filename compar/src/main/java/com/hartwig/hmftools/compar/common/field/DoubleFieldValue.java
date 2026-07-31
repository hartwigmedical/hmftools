package com.hartwig.hmftools.compar.common.field;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import java.util.List;

public class DoubleFieldValue extends FieldValue
{
    public final double Value;

    public DoubleFieldValue(final String name, final double value, final ThresholdFieldCheck check, final String formatString)
    {
        super(name, check, formatString);
        Value = value;
    }

    @Override
    public boolean checkDifference(final FieldValue other)
    {
        DoubleFieldValue otherDouble = (DoubleFieldValue)other;

        if(Value == otherDouble.Value)
            return false;

        ThresholdFieldCheck thresholdFieldCheck = (ThresholdFieldCheck)Check;

        double absDiff = abs(Value - otherDouble.Value);
        double percDiff = absDiff / max(abs(Value), abs(otherDouble.Value));

        boolean satisfiesAbsDiff = thresholdFieldCheck.AbsoluteDiff == null || absDiff > thresholdFieldCheck.AbsoluteDiff;
        boolean satisfiesRelDiff = thresholdFieldCheck.PercentageDiff == null || percDiff > thresholdFieldCheck.PercentageDiff;

        return satisfiesAbsDiff && satisfiesRelDiff;
    }

    @Override
    public String toString() { return format("%s=%s", Name, format(FormatString, Value)); }

    @Override
    public String displayValue() { return format("%s", format(FormatString, Value)); }
}
