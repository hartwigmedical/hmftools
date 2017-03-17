package com.hartwig.hmftools.patientreporter.algo;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum NotSequenceableReason {
    LOW_TUMOR_PERCENTAGE("low_tumor_percentage"),
    LOW_DNA_YIELD("low_dna_yield"),
    OTHER(Strings.EMPTY);

    @NotNull
    private final String identifier;

    NotSequenceableReason(@NotNull final String identifier) {
        this.identifier = identifier;
    }

    @NotNull
    String identifier() {
        return identifier;
    }

    @NotNull
    public static NotSequenceableReason fromIdentifier(@Nullable final String identifier) {
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
