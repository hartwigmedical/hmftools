package com.hartwig.hmftools.common.exception;

import org.jetbrains.annotations.NotNull;

public class MalformedFileException extends HartwigException {

    public MalformedFileException(@NotNull final String message) {
        super(message);
    }
}
