package com.hartwig.hmftools.datamodel.variant.msi;

import org.jetbrains.annotations.NotNull;

public enum MicrosatelliteStatus
{
    MSI("Unstable"),
    MSS("Stable"),
    UNKNOWN("Unknown");

    public static final double MSI_THRESHOLD = 4.0;

    private final String display;

    MicrosatelliteStatus(final String display)
    {
        this.display = display;
    }

    @NotNull
    public String display()
    {
        return display;
    }
}
