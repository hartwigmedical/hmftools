package com.hartwig.hmftools.common.virus;

import org.jetbrains.annotations.NotNull;

public enum VirusConstants
{
    MCV("MCV"),
    EBV("EBV"),
    HPV("HPV"),
    HBV("HBV"),
    HHV8("HHV-8");

    @NotNull
    private final String virusName;

    VirusConstants(@NotNull final String virusName)
    {
        this.virusName = virusName;
    }

    @NotNull
    public static VirusConstants fromVirusName(@NotNull String virusName)
    {
        for(final VirusConstants value : VirusConstants.values())
        {
            if(value.toString().equals(virusName))
            {
                return value;
            }
        }
        throw new IllegalStateException("Cannot resolve virus name: " + virusName);
    }

    @Override
    public String toString()
    {
        return virusName;
    }
}