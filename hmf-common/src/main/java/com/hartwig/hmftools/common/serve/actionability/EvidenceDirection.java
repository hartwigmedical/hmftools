package com.hartwig.hmftools.common.serve.actionability;

public enum EvidenceDirection {
    RESPONSIVE(true),
    PREDICTED_RESPONSIVE(false),
    RESISTANT(true),
    PREDICTED_RESISTANT(false);

    private final boolean isCertain;

    EvidenceDirection(final boolean isCertain) {
        this.isCertain = isCertain;
    }

    public boolean isCertain() {
        return isCertain;
    }
}
