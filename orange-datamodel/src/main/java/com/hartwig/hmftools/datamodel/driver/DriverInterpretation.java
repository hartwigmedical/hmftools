package com.hartwig.hmftools.datamodel.driver;

import org.jetbrains.annotations.NotNull;

public enum DriverInterpretation
{
    // must be listed from low to high for comparison
    LOW,
    MEDIUM,
    HIGH;

    @NotNull
    public static DriverInterpretation interpret(double driverLikelihood)
    {
        if(driverLikelihood > 0.8)
        {
            return HIGH;
        }
        else if(driverLikelihood > 0.2)
        {
            return MEDIUM;
        }
        else
        {
            return LOW;
        }
    }
}
