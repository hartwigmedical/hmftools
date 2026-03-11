package com.hartwig.hmftools.finding.datamodel;

import java.util.Comparator;

import jakarta.validation.constraints.NotNull;

public enum DriverInterpretation
{
    // must be listed from low to high for comparison
    UNKNOWN,
    LOW,
    MEDIUM,
    HIGH;

    public static final double DRIVER_LIKELIHOOD_LOW_THRESHOLD = 0.2;
    public static final double DRIVER_LIKELIHOOD_MEDIUM_THRESHOLD = 0.8;

    @NotNull
    public static DriverInterpretation interpret(double driverLikelihood)
    {
        if(driverLikelihood >= DRIVER_LIKELIHOOD_MEDIUM_THRESHOLD)
            return HIGH;

        if(driverLikelihood >= DRIVER_LIKELIHOOD_LOW_THRESHOLD)
            return MEDIUM;

        return LOW;
    }

    public static Comparator<DriverInterpretation> highToLow() {
        return Comparator.comparing(DriverInterpretation::ordinal).reversed();
    }
}
