package com.hartwig.hmftools.common.chord;

import org.jetbrains.annotations.NotNull;

public enum ChordStatus {
    CANNOT_BE_DETERMINED("Cannot be determined"),
    HR_PROFICIENT("Proficient"),
    HR_DEFICIENT("Deficient"),
    UNKNOWN("Unknown");

    public static final double HRD_THRESHOLD = 0.5;

    @NotNull
    private final String display;

    ChordStatus(@NotNull final String display) {
        this.display = display;
    }

    @NotNull
    public String display() {
        return display;
    }
}
