package com.hartwig.hmftools.common.exception;

import org.jetbrains.annotations.NotNull;

public class HartwigException extends Exception {

    private static final long serialVersionUID = -7416221606300503967L;

    public HartwigException(@NotNull final String message) {
        super(message);
    }
}
