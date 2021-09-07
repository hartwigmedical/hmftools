package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;

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

    @NotNull
    public static DriverInterpretation interpret(double driverLikelihood) {
        if (driverLikelihood > 0.8) {
            return HIGH;
        } else if (driverLikelihood > 0.2) {
            return MEDIUM;
        } else {
            return LOW;
        }
    }
}
