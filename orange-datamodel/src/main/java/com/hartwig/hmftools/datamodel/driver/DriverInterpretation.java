package com.hartwig.hmftools.datamodel.driver;

import java.util.Comparator;

import org.jetbrains.annotations.NotNull;

public enum DriverInterpretation
{
    // must be listed from low to high for comparison
    NONE,
    LOW,
    MEDIUM,
    HIGH;

    // TODO remove this function
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

    public static Comparator<DriverInterpretation> highToLow() {
        return Comparator.comparing(DriverInterpretation::ordinal).reversed();
    }
}
