package com.hartwig.hmftools.patientreporter.cfreport;

import org.jetbrains.annotations.NotNull;

import java.text.DecimalFormat;

public class DataUtility {

    // Data ranges and thresholds
    public static final double IMPLIED_TUMOR_PURITY_MIN = 0;
    public static final double IMPLIED_TUMOR_PURITY_MAX = 1;




    // Number formatting
    private static final DecimalFormat PERCENTAGE_FORMAT = new DecimalFormat("#'%'");



    /**
     * Remap v from [inMin, inMax] to [0, 100]
     */
    public static double mapPercentage(final double v, final double inMin, final double inMax) {
        return map(v, inMin, inMax, 0, 100);
    }

    /**
     * Remap v from [inMin, inMax] to [outMin, outMax]
     */
    public static double map(final double v, final double inMin, final double inMax, final double outMin, final double outMax) {
        return (v - inMin) * (outMax - outMin) / (inMax - inMin) + outMin;
    }

    @NotNull
    public static String formatPercentage(final double percentage) {
        return PERCENTAGE_FORMAT.format(percentage);
    }


    /**
     * All members dealing with Tumor Mutational Load
     */
    public static class MutationalLoad {

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

    }

    /**
     * All members dealing with Tumor Mutational Load
     */
    public static class MicroSatellite {

        public static final double THRESHOLD = 4D;

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

    }


}
