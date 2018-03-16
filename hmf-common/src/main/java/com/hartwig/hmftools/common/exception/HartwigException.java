package com.hartwig.hmftools.common.exception;

import org.jetbrains.annotations.NotNull;

public class HartwigException extends Exception {

    public HartwigException(@NotNull final String message) {
        super(message);
    }
}
