package com.hartwig.hmftools.finding;

import static com.hartwig.hmftools.common.purple.PurpleQCStatus.MIN_PURITY;
import static com.hartwig.hmftools.common.purple.TumorMutationalStatus.TMB_THRESHOLD;
import static com.hartwig.hmftools.common.purple.TumorMutationalStatus.TML_THRESHOLD;
import static com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus.MSI_THRESHOLD;

import com.hartwig.hmftools.datamodel.orange.ExperimentType;
import com.hartwig.hmftools.finding.datamodel.ThresholdValue;

import org.jetbrains.annotations.NotNull;

class ThresholdValueFactory
{
    // TODO: Is this defined elsewhere?
    // TODO: Make sure LOW_PURITY status is correct for targeted mode
    // This is public as it's used in OncoAct for legacy panel reports.
    public static final double TARGETED_MIN_PURITY = 0.1;
    // TODO: Is this defined elsewhere?
    private static final double HRD_THRESHOLD = 0.5;

    @NotNull
    static ThresholdValue purityValue(double value, ExperimentType experimentType)
    {
        return new ThresholdValue(value, experimentType == ExperimentType.WHOLE_GENOME ? MIN_PURITY : TARGETED_MIN_PURITY);
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
