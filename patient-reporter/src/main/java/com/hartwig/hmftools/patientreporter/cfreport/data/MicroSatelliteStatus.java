package com.hartwig.hmftools.patientreporter.cfreport.data;

import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;

import org.jetbrains.annotations.NotNull;

public class MicroSatelliteStatus {

    public static final double RANGE_MIN = 1;
    public static final double RANGE_MAX = 100;
    public static final double THRESHOLD = MicrosatelliteStatus.MSI_THRESHOLD;

    @NotNull
    public static String interpret(double microSatelliteIndelsPerMb) {
        return MicrosatelliteStatus.fromIndelsPerMb(microSatelliteIndelsPerMb).display();
    }
}
