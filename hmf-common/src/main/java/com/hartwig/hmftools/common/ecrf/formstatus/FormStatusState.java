package com.hartwig.hmftools.common.ecrf.formstatus;

import org.jetbrains.annotations.NotNull;

public enum FormStatusState {
    SAVED("saved", 0),
    SUBMITTED("submitted", 1),
    SUBMITTED_WITH_DISCREPANCIES("submitted with discrepancies", 2),
    SUBMITTED_WITH_MISSING("submitted with missing", 3),
    VERIFIED("verified", 4),
    UNKNOWN("unknown", -1);

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

    int rank() {
        return rank;
    }

    @NotNull
    public static FormStatusState best(@NotNull FormStatusState... states) {
        assert states.length >= 1;

        FormStatusState currentBest = states[0];
        for (int i = 1; i < states.length; i++) {
            if (states[i].rank > currentBest.rank) {
                currentBest = states[i];
            }
        }
        return currentBest;
    }
}
