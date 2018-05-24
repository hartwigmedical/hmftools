package com.hartwig.hmftools.patientreporter.algo;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum NotAnalysableReason {
    LOW_TUMOR_PERCENTAGE("low_tumor_percentage"),
    LOW_DNA_YIELD("low_dna_yield"),
    POST_ANALYSIS_FAIL("post_analysis_fail"),
    UNDEFINED(Strings.EMPTY);

    @NotNull
    private final String identifier;

    NotAnalysableReason(@NotNull final String identifier) {
        this.identifier = identifier;
    }

    @NotNull
    String identifier() {
        return identifier;
    }

    @NotNull
    public static NotAnalysableReason fromIdentifier(@Nullable final String identifier) {
        if (identifier == null) {
            return UNDEFINED;
        }

        for (final NotAnalysableReason reason : NotAnalysableReason.values()) {
            if (reason.identifier().equals(identifier)) {
                return reason;
            }
        }

        return UNDEFINED;
    }
}
