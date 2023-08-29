package com.hartwig.hmftools.datamodel.purple;

public enum PurpleTumorMutationalStatus {
    HIGH("High"),
    LOW("Low"),
    UNKNOWN("Unknown");

    private final String mDisplay;

    PurpleTumorMutationalStatus(final String display)
    {
        mDisplay = display;
    }

    public String display()
    {
        return mDisplay;
    }
}
