package com.hartwig.hmftools.common.exception;

import org.jetbrains.annotations.NotNull;

public class LineNotFoundException extends HealthChecksException {

    private static final long serialVersionUID = -8396650626359037492L;

    private static final String LINE_NOT_FOUND_ERROR = "File %s does not contain lines with value %s";

    public LineNotFoundException(@NotNull final String filePath, @NotNull final String filter) {
        super(String.format(LINE_NOT_FOUND_ERROR, filePath, filter));
    }
}
