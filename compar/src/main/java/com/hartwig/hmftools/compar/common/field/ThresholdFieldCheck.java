package com.hartwig.hmftools.compar.common.field;

import static java.lang.String.format;

public class ThresholdFieldCheck extends FieldCheck
{
    public final Double AbsoluteDiff;
    public final Double PercentageDiff;

    public ThresholdFieldCheck(final boolean isCompared, final Double absoluteDiff, final Double percentageDiff)
    {
        super(isCompared);
        AbsoluteDiff = absoluteDiff;
        PercentageDiff = percentageDiff;
    }

    public String toString() { return format("active(%s) thresholds(abs=%.4f, perc=%.4f)", IsCompared, AbsoluteDiff, PercentageDiff); }
}
