package com.hartwig.hmftools.protect.report.chord;

import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

public enum ChordStatus {
    HR_DEFICIENT("HR-deficient"),
    HR_PROFICIENT("HR-proficient");

    public static final double CHORD_THRESHOLD = 0.5;

    private final String status;

    ChordStatus(final String status) {
        this.status = status;
    }

    @NotNull
    public static String status() {
        return status();
    }

    @NotNull
    public static ChordStatus formChord(double chordValue) {
        return Doubles.greaterThan(chordValue, CHORD_THRESHOLD) ? HR_DEFICIENT : HR_PROFICIENT;
    }
}
