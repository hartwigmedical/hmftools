package com.hartwig.hmftools.common.io.exception;

import java.io.IOException;

import org.jetbrains.annotations.NotNull;

public class NotFolderException extends IOException {

    private static final String MESSAGE = "Path %s is not a folder";

    public NotFolderException(@NotNull final String folderName) {
        super(String.format(MESSAGE, folderName));
    }
}
