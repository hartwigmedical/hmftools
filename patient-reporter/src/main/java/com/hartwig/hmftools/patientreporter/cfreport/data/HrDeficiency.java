package com.hartwig.hmftools.patientreporter.cfreport.data;

import org.jetbrains.annotations.NotNull;

import java.text.DecimalFormat;

/**
 * All code dealing with HrDeficiencydata for presentation
 */
public class HrDeficiency {

    public static final double RANGE_MIN = 0;
    public static final double RANGE_MAX = 1;

    /**
     * Interpret HR Deficiency value. It's assumed the purity fit has been checked
     */
    @NotNull
    public static String interpretToString(final double chordHrdScore) {
        return new DecimalFormat("#.##").format(chordHrdScore);
    }

}
