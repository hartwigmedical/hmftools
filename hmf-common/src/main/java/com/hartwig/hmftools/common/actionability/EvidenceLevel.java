package com.hartwig.hmftools.common.actionability;

import org.jetbrains.annotations.NotNull;

public enum EvidenceLevel {
    LEVEL_A("A", true),
    LEVEL_B("B", true),
    LEVEL_C("C", false),
    LEVEL_D("D", false),
    LEVEL_E("E", false);

    @NotNull
    private final String readableString;

    private final boolean includeForReporting;

    EvidenceLevel(@NotNull final String readableString, final boolean includeForReporting) {
        this.readableString = readableString;
        this.includeForReporting = includeForReporting;

    }

    @NotNull
    public String readableString() {
        return readableString;
    }

    public boolean includeInReport() {
        return includeForReporting;
    }

    @NotNull
    public static EvidenceLevel fromString(@NotNull String level) {
        switch (level.toUpperCase()) {
            case "A":
                return LEVEL_A;
            case "B":
                return LEVEL_B;
            case "C":
                return LEVEL_C;
            case "D":
                return LEVEL_D;
            case "E":
                return LEVEL_E;
            default:
                    throw new IllegalArgumentException("Unrecognized evidence level item: " + level);
        }
    }
}
