package com.hartwig.hmftools.healthchecker.io;

import com.hartwig.hmftools.common.exception.HartwigException;

import org.jetbrains.annotations.NotNull;

class NotFolderException extends HartwigException {

    private static final long serialVersionUID = -8396650626359037492L;

    private static final String MESSAGE = "Path %s is not a folder";

    NotFolderException(@NotNull final String folderName) {
        super(String.format(MESSAGE, folderName));
    }
}
