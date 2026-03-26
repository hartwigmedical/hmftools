package com.hartwig.hmftools.finding;

import static com.hartwig.hmftools.common.purple.PurpleQCStatus.MIN_PURITY;
import static com.hartwig.hmftools.common.purple.MicrosatelliteStatus.MSI_THRESHOLD;
import static com.hartwig.hmftools.common.purple.TumorMutationalStatus.TMB_THRESHOLD;
import static com.hartwig.hmftools.common.purple.TumorMutationalStatus.TML_THRESHOLD;

import com.hartwig.hmftools.finding.datamodel.ThresholdValue;

import org.jetbrains.annotations.NotNull;

class ThresholdValueFactory
{
    private static final double HRD_THRESHOLD = 0.5;

    @NotNull
    static ThresholdValue purityValue(double value)
    {
        return new ThresholdValue(value, MIN_PURITY);
    }

    @NotNull
    static ThresholdValue tmlValue(double value)
    {
        return new ThresholdValue(value, TML_THRESHOLD);
    }

    @NotNull
    static ThresholdValue hrdValue(double value)
    {
        return new ThresholdValue(value, HRD_THRESHOLD);
    }

    @NotNull
    static ThresholdValue msiValue(double value)
    {
        return new ThresholdValue(value, MSI_THRESHOLD);
    }

    @NotNull
    static ThresholdValue tmbValue(double value)
    {
        return new ThresholdValue(value, TMB_THRESHOLD);
    }
}
