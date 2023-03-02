package com.hartwig.hmftools.datamodel.purple;

public enum TumorMutationalStatus
{
    HIGH("High"),
    LOW("Low"),
    UNKNOWN("Unknown");

    private final String mDisplay;

    TumorMutationalStatus(final String display)
    {
        mDisplay = display;
    }

    public String display()
    {
        return mDisplay;
    }
}
