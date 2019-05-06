package com.hartwig.hmftools.patientreporter.cfreport.data;

import org.jetbrains.annotations.NotNull;

public final class MutationalLoad {

    public static final int RANGE_MIN = 1;
    public static final int RANGE_MAX = 1000;
    public static final int THRESHOLD = 140;

    /**
     * Interpret mutational load value to be either "High" or "Low"
     */
    @NotNull
    public static String interpretToString(final int mutationalLoad, boolean hasReliablePurityFit) {
        if (!hasReliablePurityFit) {
            return DataUtil.NAString;
        } else if (mutationalLoad > THRESHOLD) {
            return "High";
        } else {
            return "Low";
        }
    }
}
