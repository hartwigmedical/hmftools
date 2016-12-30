package com.hartwig.healthchecker.checkers.checks;

public enum SummaryMetricsCheck {
    MAPPING_PF_MISMATCH_RATE(12),
    MAPPING_PF_INDEL_RATE(14),
    MAPPING_STRAND_BALANCE(19),
    MAPPING_PCT_CHIMERA(20),
    MAPPING_PCT_ADAPTER(21);

    private final int index;

    SummaryMetricsCheck(final int index) {
        this.index = index;
    }

    public int getIndex() {
        return index;
    }
}
