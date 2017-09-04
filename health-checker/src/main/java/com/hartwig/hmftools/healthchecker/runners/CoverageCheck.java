package com.hartwig.hmftools.healthchecker.runners;

import org.jetbrains.annotations.NotNull;

public enum CoverageCheck {
    COVERAGE_10X("PCT_10X", 13),
    COVERAGE_20X("PCT_20X", 15),
    COVERAGE_30X("PCT_30X", 17),
    COVERAGE_60X("PCT_60X", 20);

    @NotNull
    private final String fieldName;
    private final int columnIndex;

    CoverageCheck(@NotNull final String fieldName, final int columnIndex) {
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
