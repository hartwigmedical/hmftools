package com.hartwig.hmftools.healthchecker.exception;

import org.jetbrains.annotations.NotNull;

public class EmptyFolderException extends HealthChecksException {

    private static final long serialVersionUID = -8396650626359037492L;

    private static final String MESSAGE = "Folder %s is empty";

    public EmptyFolderException(@NotNull final String folder) {
        super(String.format(MESSAGE, folder));
    }
}
