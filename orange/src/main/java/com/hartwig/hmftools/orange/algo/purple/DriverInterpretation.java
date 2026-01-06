package com.hartwig.hmftools.orange.algo.purple;

public enum DriverInterpretation
{
    HIGH,
    MEDIUM,
    LOW;

    public static final double DRIVER_LIKELIHOOD_LOW_THRESHOLD = 0.2;
    public static final double DRIVER_LIKELIHOOD_MEDIUM_THRESHOLD = 0.8;

    public static DriverInterpretation interpret(double driverLikelihood)
    {
        if(driverLikelihood >= DRIVER_LIKELIHOOD_MEDIUM_THRESHOLD)
            return HIGH;

        if(driverLikelihood >= DRIVER_LIKELIHOOD_LOW_THRESHOLD)
            return MEDIUM;

        return LOW;
    }
}
