package com.hartwig.hmftools.common.chord;

import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

public enum ChordStatus {
    HRD("HR Deficient"),
    HRP("HR Proficient"),
    UNKNOWN("Unknown");

    public static final double HRD_THRESHOLD = 0.5;

    private final String display;

    ChordStatus(final String display) {
        this.display = display;
    }

    @NotNull
    public static ChordStatus fromHRD(double hrd) {
        return Doubles.greaterOrEqual(hrd, HRD_THRESHOLD) ? HRD : HRP;
    }

    @NotNull
    public String display() {
        return display;
    }
}
