package com.hartwig.hmftools.patientreporter.qcfail;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum QCFailReason {
    LOW_DNA_YIELD("low_dna_yield", QCFailType.LOW_QUALITY_BIOPSY, false),
    POST_ANALYSIS_FAIL("post_analysis_fail", QCFailType.LOW_QUALITY_BIOPSY, true),
    SHALLOW_SEQ_LOW_PURITY("shallow_seq_low_purity", QCFailType.LOW_QUALITY_BIOPSY, false),
    INSUFFICIENT_TISSUE("insufficient_tissue_delivered", QCFailType.LOW_QUALITY_BIOPSY, false),
    BELOW_DETECTION_THRESHOLD("below_detection_threshold", QCFailType.LOW_QUALITY_BIOPSY, true),
    LAB_FAILURE("lab_failure", QCFailType.TECHNICAL_FAILURE, false),
    UNDEFINED(Strings.EMPTY, QCFailType.UNDEFINED, false);

    @NotNull
    private final String identifier;
    @NotNull
    private final QCFailType type;
    private final boolean fullWgsDataAvailable;

    QCFailReason(@NotNull final String identifier, @NotNull final QCFailType type, final boolean fullWgsDataAvailable) {
        this.identifier = identifier;
        this.type = type;
        this.fullWgsDataAvailable = fullWgsDataAvailable;
    }

    @NotNull
    public String identifier() {
        return identifier;
    }

    @NotNull
    public QCFailType type() {
        return type;
    }

    public boolean isFullWgsDataAvailable() {
        return fullWgsDataAvailable;
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
