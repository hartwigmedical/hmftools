package com.hartwig.hmftools.protect.report.chord;

import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

public enum ChordStatus {
    HR_DEFICIENT("HR-deficient"),
    HR_PROFICIENT("HR-proficient");

    public static final double CHORD_THRESHOLD = 0.5;

    private final String display;

    ChordStatus(final String display) {
        this.display = display;
    }

    @NotNull
    public static String display() {
        return display();
    }

    @NotNull
    public static ChordStatus formChord(double chordValue) {
        return Doubles.greaterThan(0.5, CHORD_THRESHOLD) ? HR_DEFICIENT : HR_PROFICIENT;
    }
}
