package com.hartwig.hmftools.common.actionability;

import org.jetbrains.annotations.NotNull;

public enum EvidenceLevel {
    LEVEL_A("A", true),
    LEVEL_B("B", true),
    LEVEL_C("C", false),
    LEVEL_D("D", false),
    LEVEL_E("E", false);

    @NotNull
    private final String levelEvidenceItem;

    private final boolean isReportedEvidenceItemLevel;

    EvidenceLevel(@NotNull final String levelEvidenceItem, final boolean isReportedEvidenceItemLevel) {
        this.levelEvidenceItem = levelEvidenceItem;
        this.isReportedEvidenceItemLevel = isReportedEvidenceItemLevel;

    }

    @NotNull
    public String levelEvidenceItem() {
        return levelEvidenceItem;
    }

    public boolean isReportedEvidenceItemLevel() {
        return isReportedEvidenceItemLevel;
    }

    @NotNull
    public static EvidenceLevel fromString(@NotNull String levelItem) {
        switch (levelItem.toUpperCase()) {
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
                    throw new IllegalArgumentException("Unrecognized evidence level item" + levelItem);
        }
    }
}
