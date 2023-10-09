package com.hartwig.hmftools.datamodel.purple;

import org.jetbrains.annotations.NotNull;

public enum PurpleMicrosatelliteStatus
{
    MSI("Unstable"),
    MSS("Stable"),
    UNKNOWN("Unknown");

    @NotNull
    private final String display;

    PurpleMicrosatelliteStatus(@NotNull final String display)
    {
        this.display = display;
    }

    @NotNull
    public String display()
    {
        return display;
    }
}
