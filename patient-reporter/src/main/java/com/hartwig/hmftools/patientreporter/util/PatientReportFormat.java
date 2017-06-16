package com.hartwig.hmftools.patientreporter.util;

import org.jetbrains.annotations.Nullable;

public enum PatientReportFormat {
    ;

    public static String formatPercent(final @Nullable Double percentage) {
        return formatPercent(percentage == null ? Double.NaN : percentage);
    }

    public static String formatPercent(final double percentage) {
        return Long.toString(Math.round(percentage * 100D)) + "%";
    }
}
