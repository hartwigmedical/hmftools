package com.hartwig.hmftools.protect.conclusion;

import org.jetbrains.annotations.NotNull;

public enum AberrationGenSummary {

    HIGH_MTL("High TMB"),
    HR_DEFICIENT("HR-deficient"),
    MSI("MSI"),
    LOW_PURITY("LOW PURITY");

    @NotNull
    private final String readableString;

    AberrationGenSummary(@NotNull final String readableString) {
        this.readableString = readableString;

    }

    @NotNull
    public String readableString() {
        return readableString;
    }

    @NotNull
    public static AberrationGenSummary fromString(@NotNull String gene) {
        switch (gene) {
            case "HR-deficient":
                return HR_DEFICIENT;
            case "High TMB":
                return HIGH_MTL;
            case "MSI":
                return MSI;
            case "LOW PURITY":
                return LOW_PURITY;
            default:
                throw new IllegalArgumentException("Unrecognized gene event " + gene);
        }
    }
}
