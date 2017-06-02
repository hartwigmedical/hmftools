package com.hartwig.hmftools.patientreporter.util;

public enum PatientReportFormat {
    ;

    public static String formatPercent(Double percentage) {
        return formatPercent(percentage == null ? Double.NaN : percentage);
    }

    public static String formatPercent(double percentage) {
        return Long.toString(Math.round(percentage * 100D)) + "%";
    }
}
