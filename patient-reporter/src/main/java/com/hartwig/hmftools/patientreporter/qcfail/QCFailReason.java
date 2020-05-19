package com.hartwig.hmftools.patientreporter.qcfail;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum QCFailReason {
    LOW_DNA_YIELD("low_dna_yield"),
    POST_ANALYSIS_FAIL("post_analysis_fail"),
    SHALLOW_SEQ_LOW_PURITY("shallow_seq_low_purity"),
    INSUFFICIENT_TISSUE("insufficient_tissue_delivered"),
    BELOW_DETECTION_THRESHOLD("below_detection_threshold"),
    LAB_FAILURE("lab_failure"),
    UNDEFINED(Strings.EMPTY);

    @NotNull
    private final String identifier;

    QCFailReason(@NotNull final String identifier) {
        this.identifier = identifier;
    }

    @NotNull
    public String identifier() {
        return identifier;
    }

    @NotNull
    public static QCFailReason fromIdentifier(@Nullable final String identifier) {
        if (identifier == null) {
            return UNDEFINED;
        }

        for (QCFailReason reason : QCFailReason.values()) {
            if (reason.identifier().equals(identifier)) {
                return reason;
            }
        }

        return UNDEFINED;
    }
}
