package com.hartwig.hmftools.patientreporter.qcfail;

import com.hartwig.hmftools.patientreporter.PatientReporterApplication;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum QCFailReason {
    LOW_TUMOR_PERCENTAGE("low_tumor_percentage", "HMF Low Tumor Percentage Report v" + PatientReporterApplication.VERSION),
    LOW_DNA_YIELD("low_dna_yield", "HMF Low DNA Yield Report v" + PatientReporterApplication.VERSION),
    POST_ANALYSIS_FAIL("post_analysis_fail", "HMF Failed Post DNA Isolation Report v" + PatientReporterApplication.VERSION),
    SHALLOW_SEQ_LOW_PURITY("shallow_seq_low_purity",
            "HMF Low Molecular Tumor Percentage Report v" + "" + PatientReporterApplication.VERSION),
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
    String identifier() {
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
