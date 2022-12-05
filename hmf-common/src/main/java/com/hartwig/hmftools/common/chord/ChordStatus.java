package com.hartwig.hmftools.common.chord;

public enum ChordStatus
{
    CANNOT_BE_DETERMINED("Cannot be determined"),
    HR_PROFICIENT("Proficient"),
    HR_DEFICIENT("Deficient"),
    UNKNOWN("Unknown");

    private final String display;

    ChordStatus(final String display)
    {
        this.display = display;
    }

    public String display()
    {
        return display;
    }
}
