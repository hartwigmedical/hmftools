package com.hartwig.hmftools.finding;

import com.hartwig.hmftools.finding.datamodel.ThresholdValue;

import org.jetbrains.annotations.NotNull;

class ThresholdValueFactory
{
    private static final double TML_RANGE_MIN = 1;
    private static final double TML_RANGE_MAX = 1000;
    // TODO: Lookup proper threshold constant
    private static final double TML_THRESHOLD = 10;
    private static final double HRD_RANGE_MIN = 0;
    private static final double HRD_RANGE_MAX = 1;
    // TODO: Lookup proper threshold constant
    private static final double HRD_THRESHOLD = 0.5;
    private static final double MSS_RANGE_MIN = 1;
    private static final double MSS_RANGE_MAX = 100;
    // TODO: Lookup proper threshold constant
    private static final double MSS_THRESHOLD = 4.0;
    private static final double TMB_RANGE_MIN = 1;
    private static final double TMB_RANGE_MAX = 120;
    // TODO: Lookup proper threshold constant
    private static final double TMB_THRESHOLD = 10;

    @NotNull
    static ThresholdValue tmlValue(double value)
    {
        return new ThresholdValue(value, TML_RANGE_MIN, TML_RANGE_MAX, TML_THRESHOLD);
    }

    @NotNull
    static ThresholdValue hrdValue(double value)
    {
        return new ThresholdValue(value, HRD_RANGE_MIN, HRD_RANGE_MAX, HRD_THRESHOLD);
    }

    @NotNull
    static ThresholdValue mssValue(double value)
    {
        return new ThresholdValue(value, MSS_RANGE_MIN, MSS_RANGE_MAX, MSS_THRESHOLD);
    }

    @NotNull
    static ThresholdValue tmbValue(double value)
    {
        return new ThresholdValue(value, TMB_RANGE_MIN, TMB_RANGE_MAX, TMB_THRESHOLD);
    }
}
