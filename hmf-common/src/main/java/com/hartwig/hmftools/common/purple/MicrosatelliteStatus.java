package com.hartwig.hmftools.common.purple;

import com.hartwig.hmftools.common.utils.Doubles;

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
    public static MicrosatelliteStatus fromIndelsPerMb(double microsatelliteIndelsPerMb)
    {
        return Doubles.greaterOrEqual(microsatelliteIndelsPerMb, MSI_THRESHOLD) ? MSI : MSS;
    }

    @NotNull
    public String display()
    {
        return display;
    }
}
