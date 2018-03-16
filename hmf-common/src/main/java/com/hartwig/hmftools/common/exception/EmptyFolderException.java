package com.hartwig.hmftools.common.exception;

import org.jetbrains.annotations.NotNull;

public class EmptyFolderException extends HartwigException {

    private static final String MESSAGE = "Folder %s is empty";

    public EmptyFolderException(@NotNull final String folder) {
        super(String.format(MESSAGE, folder));
    }
}
