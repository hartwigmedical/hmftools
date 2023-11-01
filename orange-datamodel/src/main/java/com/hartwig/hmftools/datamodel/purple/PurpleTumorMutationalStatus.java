package com.hartwig.hmftools.datamodel.purple;

import org.jetbrains.annotations.NotNull;

public enum PurpleTumorMutationalStatus
{
    HIGH("High"),
    LOW("Low"),
    UNKNOWN("Unknown");

    @NotNull
    private final String display;

    PurpleTumorMutationalStatus(@NotNull final String display)
    {
        this.display = display;
    }

    @NotNull
    public String display()
    {
        return display;
    }
}
