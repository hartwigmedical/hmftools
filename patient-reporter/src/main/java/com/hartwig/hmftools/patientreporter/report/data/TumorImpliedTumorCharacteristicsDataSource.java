package com.hartwig.hmftools.patientreporter.report.data;

import com.hartwig.hmftools.patientreporter.report.components.MicrosatelliteSection;

import org.jetbrains.annotations.NotNull;

public class TumorImpliedTumorCharacteristicsDataSource {

    private TumorImpliedTumorCharacteristicsDataSource() {
    }

    @NotNull
    public static String interpretMutationalLoad(int mutationalLoad) {
        if (mutationalLoad > MicrosatelliteSection.MSI_THRESHOLD) {
            return "High";
        } else {
            return "Low";
        }
    }

    @NotNull
    public static String interpretMSI(double microsatelliteIndicator) {
        if (microsatelliteIndicator > MicrosatelliteSection.MSI_THRESHOLD) {
            return "Unstable";
        } else {
            return "Stable";
        }
    }


}
