package com.hartwig.hmftools.patientreporter.copynumber;

import org.jetbrains.annotations.NotNull;

public enum CopyNumberReportType {
    GAIN("copy-gain"),
    LOSS("copy-loss"),
    NEUTRAL("none");

    private final String description;

    CopyNumberReportType(final String description) {
        this.description = description;
    }

    public String description() {
        return description;
    }

    @NotNull
    public static CopyNumberReportType resolveType(final int copyNumber) {
        if (copyNumber > 2) {
            return GAIN;
        } else if (copyNumber < 2) {
            return LOSS;
        } else {
            return NEUTRAL;
        }
    }
}
