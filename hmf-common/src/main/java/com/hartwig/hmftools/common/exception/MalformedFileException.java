package com.hartwig.hmftools.common.exception;

import java.io.IOException;

import org.jetbrains.annotations.NotNull;

public class MalformedFileException extends IOException {

    public MalformedFileException(@NotNull final String message) {
        super(message);
    }
}
