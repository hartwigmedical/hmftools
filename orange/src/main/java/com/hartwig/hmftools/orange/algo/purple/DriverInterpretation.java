package com.hartwig.hmftools.orange.algo.purple;

import org.jetbrains.annotations.NotNull;

public enum DriverInterpretation {
    HIGH,
    MEDIUM,
    LOW;

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
