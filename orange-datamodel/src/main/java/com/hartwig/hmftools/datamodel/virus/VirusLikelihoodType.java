package com.hartwig.hmftools.datamodel.virus;

import org.jetbrains.annotations.NotNull;

public enum VirusLikelihoodType {
    HIGH("High"),
    LOW("Low"),
    UNKNOWN("Unknown"),;

    private final String mDisplay;

    VirusLikelihoodType(@NotNull final String display) {
        mDisplay = display;
    }

    public String display() {
        return mDisplay;
    }
}
