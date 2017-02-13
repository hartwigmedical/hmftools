package com.hartwig.hmftools.patientreporter;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum NotSequenceableReason {
    LOW_TUMOR_PERCENTAGE("low_tumor_percentage", "Tumor Percentage <30% in biopsy"),
    LOW_DNA_YIELD("low_dna_yield", "Not enough DNA available from biopsy"),
    OTHER(Strings.EMPTY, Strings.EMPTY);

    @NotNull
    private final String identifier;
    @NotNull
    private final String message;

    NotSequenceableReason(@NotNull final String identifier, @NotNull final String message) {
        this.identifier = identifier;
        this.message = message;
    }

    @NotNull
    String identifier() {
        return identifier;
    }

    @NotNull
    public String message() {
        return message;
    }

    @NotNull
    static NotSequenceableReason fromIdentifier(@Nullable final String identifier) {
        if (identifier == null) {
            return OTHER;
        }

        for (final NotSequenceableReason reason : NotSequenceableReason.values()) {
            if (reason.identifier().equals(identifier)) {
                return reason;
            }
        }

        return OTHER;
    }
}
