package com.hartwig.hmftools.common.serve.actionability;

public enum EvidenceDirection {
    RESPONSIVE(true, true),
    PREDICTED_RESPONSIVE(true, false),
    RESISTANT(true, true),
    PREDICTED_RESISTANT(true, false);

    private final boolean isResponsive;
    private final boolean isCertain;

    EvidenceDirection(final boolean isResponsive, final boolean isCertain) {
        this.isResponsive = isResponsive;
        this.isCertain = isCertain;
    }

    public boolean isResponsive() {
        return isResponsive;
    }

    public boolean isCertain() {
        return isCertain;
    }
}
