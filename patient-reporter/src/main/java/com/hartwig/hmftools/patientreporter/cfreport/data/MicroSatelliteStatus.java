package com.hartwig.hmftools.patientreporter.cfreport.data;

import org.jetbrains.annotations.NotNull;

public enum MicroSatelliteStatus {
    UNSTABLE("Unstable"),
    STABLE("Stable");

    public static final double RANGE_MIN = 1;
    public static final double RANGE_MAX = 100;
    public static final double THRESHOLD = 4;

    @NotNull
    private final String text;

    MicroSatelliteStatus(@NotNull final String text) {
        this.text = text;
    }

    @NotNull
    public String text() {
        return text;
    }

    @NotNull
    public static MicroSatelliteStatus interpret(double microSatelliteIndelsPerMb) {
        if (microSatelliteIndelsPerMb > THRESHOLD) {
            return UNSTABLE;
        } else {
            return STABLE;
        }
    }
}
