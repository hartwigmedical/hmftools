package com.hartwig.hmftools.common.chord;

import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;

import org.jetbrains.annotations.NotNull;

public enum ChordStatus {
    CANNOT_BE_DETERMINED("Cannot be determined"),
    HR_PROFICIENT("Proficient"),
    HR_DEFICIENT("Deficient"),
    UNKNWON("Unkown");

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
