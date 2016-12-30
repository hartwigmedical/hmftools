package com.hartwig.hmftools.healthchecker.runners.checks;

import org.jetbrains.annotations.NotNull;

public enum InsertSizeMetricsCheck {
    MAPPING_MEDIAN_INSERT_SIZE("MEDIAN_INSERT_SIZE", 0),
    MAPPING_WIDTH_OF_70_PERCENT("WIDTH_OF_70_PERCENT", 14);

    @NotNull
    private final String fieldName;
    private final int columnIndex;

    InsertSizeMetricsCheck(@NotNull final String fieldName, final int columnIndex) {
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
