package com.hartwig.hmftools.patientreporter.cfreport.data;

import org.jetbrains.annotations.NotNull;

/**
 * All code dealing with Mutational Load for presentation
 */
public class MutationalLoad {

    public static final int RANGE_MIN = 1;
    public static final int RANGE_MAX = 1000;
    public static final int THRESHOLD = 140;

    /**
     * Interpret mutational load value to be either "High" or "Low". It's assumed the purity
     * fit has been checked
     */
    @NotNull
    public static String interpretToString(final int mutationalLoad) {
        if (mutationalLoad > THRESHOLD) {
            return "High";
        } else {
            return "Low";
        }
    }

    public static double scale(double value) {
        return Math.log10(value);
    }

}
