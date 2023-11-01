package com.hartwig.hmftools.datamodel.chord;

import org.jetbrains.annotations.NotNull;

public enum ChordStatus
{
    CANNOT_BE_DETERMINED("Cannot be determined"),
    HR_PROFICIENT("Proficient"),
    HR_DEFICIENT("Deficient"),
    UNKNOWN("Unknown");

    @NotNull
    private final String display;

    ChordStatus(@NotNull final String display)
    {
        this.display = display;
    }

    @NotNull
    public String display()
    {
        return display;
    }
}
