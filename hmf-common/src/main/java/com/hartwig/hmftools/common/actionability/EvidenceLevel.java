package com.hartwig.hmftools.common.actionability;

import org.jetbrains.annotations.NotNull;

public enum EvidenceLevel {
    LEVELA("A", true),
    LEVELB("B", true),
    LEVELC("C", false),
    LEVELD("D", false),
    LEVELE("E", false);

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
                return LEVELA;
            case "B":
                return LEVELB;
            case "C":
                return LEVELC;
            case "D":
                return LEVELD;
            case "E":
                return LEVELE;
            default:
                    throw new IllegalArgumentException("Unrecognized evidence level item" + levelItem);
        }
    }
}
