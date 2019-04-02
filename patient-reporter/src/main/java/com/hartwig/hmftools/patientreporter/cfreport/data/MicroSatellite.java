package com.hartwig.hmftools.patientreporter.cfreport.data;

import org.jetbrains.annotations.NotNull;

/**
 * All code dealing with Micro Satellite Stability for presentation
 */
public class MicroSatellite {

    public static final double RANGE_MIN = 1E-2;
    public static final double RANGE_MAX = 100;
    public static final double THRESHOLD = 4;

    /**
     * Interpret micro satellite indel value to be either "Stable" or "Instable". It's assumed the purity
     * fit has been checked
     */
    @NotNull
    public static String interpretToString(double microSatelliteIndelsPerMb) {
        if (microSatelliteIndelsPerMb > THRESHOLD) {
            return "Unstable";
        } else {
            return "Stable";
        }
    }

    public static double scale(double value) {
        return Math.log10(value);
    }

}
