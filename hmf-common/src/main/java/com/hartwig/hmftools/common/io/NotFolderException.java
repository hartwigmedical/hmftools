package com.hartwig.hmftools.common.io;

import java.io.IOException;

import org.jetbrains.annotations.NotNull;

class NotFolderException extends IOException {

    private static final String MESSAGE = "Path %s is not a folder";

    NotFolderException(@NotNull final String folderName) {
        super(String.format(MESSAGE, folderName));
    }
}
