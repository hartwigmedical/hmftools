package com.hartwig.hmftools.patientreporter.cfreport.data;

import org.jetbrains.annotations.NotNull;

import java.text.DecimalFormat;

public final class HrDeficiency {

    public static final double RANGE_MIN = 0;
    public static final double RANGE_MAX = 1;

    /**
     * Interpret HR Deficiency value
     */
    @NotNull
    public static String interpretToString(final double chordHrdScore, boolean hasReliablePurityFit) {
        if (!hasReliablePurityFit) {
            return DataUtil.NAString;
        } else {
            return new DecimalFormat("#.##").format(chordHrdScore);
        }
    }
}
