package com.hartwig.hmftools.common.chord;

import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;

import org.jetbrains.annotations.NotNull;

public enum ChordStatus {
    CANNOT_BE_DETERMINED("cannot_be_determined", "Cannot be determined"),
    HR_PROFICIENT("HR_proficient", "Proficient"),
    HR_DEFICIENT("HR_deficient", "Deficient"),
    UNKNOWN("unkown", "UnkNown");

    public static final double HRD_THRESHOLD = 0.5;
    @NotNull
    private final String status;
    @NotNull
    private final String display;

    ChordStatus(@NotNull final String status, @NotNull final String display) {
        this.status = status;
        this.display = display;
    }

    @NotNull
    public String status() {
        return status;
    }
    @NotNull
    public String display() {
        return display;
    }
}
