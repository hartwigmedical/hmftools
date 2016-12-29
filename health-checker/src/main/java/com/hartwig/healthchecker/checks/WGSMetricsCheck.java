package com.hartwig.healthchecker.checks;

import org.jetbrains.annotations.NotNull;

enum WGSMetricsCheck {
    COVERAGE_MEAN("MEAN_COVERAGE", 1),
    COVERAGE_MEDIAN("MEDIAN_COVERAGE", 3),
    COVERAGE_SD("SD_COVERAGE", 2),
    COVERAGE_10X("PCT_10X", 13),
    COVERAGE_20X("PCT_20X", 15),
    COVERAGE_30X("PCT_30X", 17),
    COVERAGE_60X("PCT_60X", 20),
    COVERAGE_PCT_EXC_MAPQ("PCT_EXC_MAPQ", 5),
    COVERAGE_PCT_EXC_DUPE("PCT_EXC_DUPE", 6),
    COVERAGE_PCT_EXC_UNPAIRED("PCT_EXC_UNPAIRED", 7),
    COVERAGE_PCT_EXC_BASEQ("PCT_EXC_BASEQ", 8),
    COVERAGE_PCT_EXC_OVERLAP("PCT_EXC_OVERLAP", 9),
    COVERAGE_PCT_EXC_TOTAL("PCT_EXC_TOTAL", 11),;

    @NotNull
    private final String fieldName;
    private final int columnIndex;

    WGSMetricsCheck(@NotNull final String fieldName, final int columnIndex) {
        this.fieldName = fieldName;
        this.columnIndex = columnIndex;
    }

    @NotNull
    public String getFieldName() {
        return fieldName;
    }

    public int getColumnIndex() {
        return columnIndex;
    }
}
