package com.hartwig.hmftools.patientreporter.algo;

import com.hartwig.hmftools.patientreporter.PatientReporterApplication;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum NotAnalysableReason {
    LOW_TUMOR_PERCENTAGE("low_tumor_percentage", "HMF Low Tumor Percentage Report v" + PatientReporterApplication.VERSION),
    LOW_DNA_YIELD("low_dna_yield", "HMF Low DNA Yield Report v" + PatientReporterApplication.VERSION),
    POST_ANALYSIS_FAIL("post_analysis_fail", "HMF Failed Post DNA Isolation Report v" + PatientReporterApplication.VERSION),
    UNDEFINED(Strings.EMPTY, Strings.EMPTY);

    @NotNull
    private final String identifier;
    @NotNull
    private final String title;

    NotAnalysableReason(@NotNull final String identifier, @NotNull final String title) {
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
