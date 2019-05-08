package com.hartwig.hmftools.patientreporter.cfreport.data;

import org.jetbrains.annotations.NotNull;

import java.text.DecimalFormat;

public final class HrDeficiency {

    public static final double RANGE_MIN = 0;
    public static final double RANGE_MAX = 1;

    private HrDeficiency() {
    }

    @NotNull
    public static String interpretToString(final double chordHrdScore, boolean hasReliablePurityFit) {
        if (!hasReliablePurityFit) {
            return DataUtil.NA_STRING;
        } else {
            return new DecimalFormat("#.##").format(chordHrdScore);
        }
    }
}
