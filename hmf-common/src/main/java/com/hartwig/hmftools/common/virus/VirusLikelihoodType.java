package com.hartwig.hmftools.common.virus;

import org.jetbrains.annotations.NotNull;

public enum VirusLikelihoodType {
    HIGH("HIGH"),
    LOW("LOW"),
    UNKNOWN("UNKNOWN"),;

    private final String mDisplay;

    VirusLikelihoodType(@NotNull final String display) {
        mDisplay = display;
    }

    public String display() {
        return mDisplay;
    }
}
