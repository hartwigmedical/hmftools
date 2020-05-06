package com.hartwig.hmftools.patientreporter.qcfail;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum QCFailReason {
    LOW_DNA_YIELD("low_dna_yield", "Low DNA Yield Report"),
    POST_ANALYSIS_FAIL("post_analysis_fail", "Failed Post DNA Isolation Report"),
    SHALLOW_SEQ_LOW_PURITY("shallow_seq_low_purity", "Low Molecular Tumor Percentage Report"),
    INSUFFICIENT_TISSUE("insufficient_tissue_delivered", "Insufficient Tissue Report"),
    BELOW_DETECTION_THRESHOLD("below_detection_threshold", "No Observed Variants Report"),
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
