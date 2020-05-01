package com.hartwig.hmftools.patientreporter.cfreport.data;

import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;

import org.jetbrains.annotations.NotNull;

public final class HrDeficiency {

    public static final double RANGE_MIN = 0;
    public static final double RANGE_MAX = 1;

    private HrDeficiency() {
    }

    @NotNull
    public static String interpretToString(final double chordHrdScore) {
        return ReportResources.decimalFormat("#.##").format(chordHrdScore);
    }
}
