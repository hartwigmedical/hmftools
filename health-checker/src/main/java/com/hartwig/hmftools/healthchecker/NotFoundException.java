package com.hartwig.hmftools.healthchecker;

import com.hartwig.hmftools.common.exception.HartwigException;

import org.jetbrains.annotations.NotNull;

class NotFoundException extends HartwigException {

    private static final long serialVersionUID = 6625914789247983088L;

    NotFoundException(@NotNull final String message) {
        super(message);
    }
}
