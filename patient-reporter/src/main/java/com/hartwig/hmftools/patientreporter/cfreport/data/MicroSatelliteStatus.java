package com.hartwig.hmftools.patientreporter.cfreport.data;

import org.jetbrains.annotations.NotNull;

public final class MicroSatelliteStatus {

    public static final double RANGE_MIN = 1;
    public static final double RANGE_MAX = 100;
    public static final double THRESHOLD = 4;

    /**
     * Interpret micro satellite indel value to be either "Stable" or "Instable"
     */
    @NotNull
    public static String interpretToString(double microSatelliteIndelsPerMb) {
        if (microSatelliteIndelsPerMb > THRESHOLD) {
            return "Unstable";
        } else {
            return "Stable";
        }
    }
}
