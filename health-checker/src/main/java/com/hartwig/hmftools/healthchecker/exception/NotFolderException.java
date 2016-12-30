package com.hartwig.hmftools.healthchecker.exception;

import org.jetbrains.annotations.NotNull;

public class NotFolderException extends HealthChecksException {

    private static final long serialVersionUID = -8396650626359037492L;

    private static final String MESSAGE = "Path %s is not a folder";

    public NotFolderException(@NotNull final String folderName) {
        super(String.format(MESSAGE, folderName));
    }
}
