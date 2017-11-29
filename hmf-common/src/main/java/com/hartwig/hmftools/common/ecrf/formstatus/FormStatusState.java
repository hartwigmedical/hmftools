package com.hartwig.hmftools.common.ecrf.formstatus;

import org.jetbrains.annotations.NotNull;

public enum FormStatusState {
    SAVED("saved"),
    SUBMITTED("submitted"),
    SUBMITTED_WITH_DISCREPANCIES("submitted with discrepancies"),
    SUBMITTED_WITH_MISSING("submitted with missing"),
    VERIFIED("verified"),
    UNKNOWN("unknown");

    @NotNull
    private final String stateString;

    FormStatusState(@NotNull final String stateString) {
        this.stateString = stateString;
    }

    @NotNull
    public String stateString() {
        return stateString;
    }
}
