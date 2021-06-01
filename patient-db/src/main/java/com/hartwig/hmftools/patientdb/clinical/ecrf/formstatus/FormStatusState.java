package com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus;

import org.jetbrains.annotations.NotNull;

public enum FormStatusState {
    SAVED("saved", 0),
    SUBMITTED("submitted", 1),
    SUBMITTED_WITH_DISCREPANCIES("submitted with discrepancies", 2),
    SUBMITTED_WITH_MISSING("submitted with missing", 3),
    VERIFIED("verified", 4),
    UNDEFINED("undefined", -1);

    @NotNull
    private final String stateString;
    private final int rank;

    FormStatusState(@NotNull final String stateString, final int rank) {
        this.stateString = stateString;
        this.rank = rank;
    }

    @NotNull
    public String stateString() {
        return stateString;
    }

    @NotNull
    public static FormStatusState best(@NotNull FormStatusState state1, @NotNull FormStatusState state2) {
        return state1.rank > state2.rank ? state1 : state2;
    }
}
