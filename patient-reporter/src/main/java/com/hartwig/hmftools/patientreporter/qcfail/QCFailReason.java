package com.hartwig.hmftools.patientreporter.qcfail;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum QCFailReason {
    LOW_DNA_YIELD("low_dna_yield", "Failed Sample Report"),
    POST_ANALYSIS_FAIL("post_analysis_fail", "Failed Sample Report"),
    SHALLOW_SEQ_LOW_PURITY("shallow_seq_low_purity", "Failed Sample Report"),
    INSUFFICIENT_TISSUE("insufficient_tissue_delivered", "Failed Sample Report"),
    BELOW_DETECTION_THRESHOLD("below_detection_threshold", "Failed Sample Report"),
    LAB_FAILURE("lab_failure", "Failed Sample Report"),
    UNDEFINED(Strings.EMPTY, Strings.EMPTY);

    @NotNull
    private final String identifier;
    @NotNull
    private final String title;

    QCFailReason(@NotNull final String identifier, @NotNull final String title) {
        this.identifier = identifier;
        this.title = title;
    }

    @NotNull
    public String identifier() {
        return identifier;
    }

    @NotNull
    public String title() {
        return title;
    }

    @NotNull
    public static QCFailReason fromIdentifier(@Nullable final String identifier) {
        if (identifier == null) {
            return UNDEFINED;
        }

        for (final QCFailReason reason : QCFailReason.values()) {
            if (reason.identifier().equals(identifier)) {
                return reason;
            }
        }

        return UNDEFINED;
    }
}
