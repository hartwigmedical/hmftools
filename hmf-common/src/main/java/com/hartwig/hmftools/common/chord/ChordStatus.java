package com.hartwig.hmftools.common.chord;

import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;

import org.jetbrains.annotations.NotNull;

public enum ChordStatus {
    HRD("Deficient"),
    HRP("Proficient"),
    UNKNOWN("Unknown");

    public static final double HRD_THRESHOLD = 0.5;

    @NotNull
    private final String display;

    ChordStatus(@NotNull final String display) {
        this.display = display;
    }

    @NotNull
    public static ChordStatus fromHRD(double hrd) {
        return Doubles.greaterOrEqual(hrd, HRD_THRESHOLD) ? HRD : HRP;
    }

    public static boolean reliableHRD(@NotNull MicrosatelliteStatus microsatelliteStatus) {
        return microsatelliteStatus == MicrosatelliteStatus.MSS;
    }

    @NotNull
    public String display() {
        return display;
    }
}
