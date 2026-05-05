package com.hartwig.hmftools.finding;

import static com.hartwig.hmftools.common.purple.PurpleQCStatus.MIN_PURITY;
import static com.hartwig.hmftools.common.purple.TumorMutationalStatus.TMB_THRESHOLD;
import static com.hartwig.hmftools.common.purple.TumorMutationalStatus.TML_THRESHOLD;
import static com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus.MSI_THRESHOLD;

import com.hartwig.hmftools.datamodel.orange.ExperimentType;
import com.hartwig.hmftools.finding.datamodel.ThresholdValue;

class ThresholdValueFactory
{
    // TODO: Is this defined elsewhere?
    private static final double HRD_THRESHOLD = 0.5;

    static ThresholdValue tmlValue(double value)
    {
        return new ThresholdValue(value, TML_THRESHOLD);
    }

    static ThresholdValue hrdValue(double value)
    {
        return new ThresholdValue(value, HRD_THRESHOLD);
    }

    static ThresholdValue msiValue(double value)
    {
        return new ThresholdValue(value, MSI_THRESHOLD);
    }

    static ThresholdValue tmbValue(double value)
    {
        return new ThresholdValue(value, TMB_THRESHOLD);
    }
}
