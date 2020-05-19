package com.hartwig.hmftools.patientreporter.qcfail;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum QCFailReason {
    LOW_DNA_YIELD("low_dna_yield", QCFailType.LOW_QUALITY_BIOPSY),
    POST_ANALYSIS_FAIL("post_analysis_fail", QCFailType.LOW_QUALITY_BIOPSY),
    SHALLOW_SEQ_LOW_PURITY("shallow_seq_low_purity", QCFailType.LOW_QUALITY_BIOPSY),
    INSUFFICIENT_TISSUE("insufficient_tissue_delivered", QCFailType.LOW_QUALITY_BIOPSY),
    BELOW_DETECTION_THRESHOLD("below_detection_threshold", QCFailType.LOW_QUALITY_BIOPSY),
    LAB_FAILURE("lab_failure", QCFailType.TECHNICAL_FAILURE),
    UNDEFINED(Strings.EMPTY, QCFailType.UNDEFINED);

    @NotNull
    private final String identifier;
    @NotNull
    private final QCFailType type;

    QCFailReason(@NotNull final String identifier, @NotNull final QCFailType type) {
        this.identifier = identifier;
        this.type = type;
    }

    @NotNull
    public String identifier() {
        return identifier;
    }

    @NotNull
    public QCFailType type() {
        return type;
    }

    @NotNull
    public static QCFailReason fromIdentifier(@Nullable String identifier) {
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
