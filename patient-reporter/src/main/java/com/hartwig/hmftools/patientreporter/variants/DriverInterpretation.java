package com.hartwig.hmftools.patientreporter.variants;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum DriverInterpretation {
    HIGH("High"),
    MEDIUM("Medium"),
    LOW("Low");

    @NotNull
    private final String text;

    DriverInterpretation(@NotNull final String text) {
        this.text = text;
    }

    @NotNull
    public String text() {
        return text;
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
