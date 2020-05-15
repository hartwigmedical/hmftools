package com.hartwig.hmftools.patientreporter.variants;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum DriverInterpretation {
    HIGH("High"),
    MEDIUM("Medium"),
    LOW("Low");

    @NotNull
    private final String display;

    DriverInterpretation(@NotNull final String display) {
        this.display = display;
    }

    @NotNull
    public String display() {
        return display;
    }

    @Nullable
    public static DriverInterpretation interpret(@Nullable Double driverLikelihood) {
        if (driverLikelihood == null) {
            return null;
        }

        if (driverLikelihood > 0.8) {
            return HIGH;
        } else if (driverLikelihood > 0.2) {
            return MEDIUM;
        } else {
            return LOW;
        }
    }
}
