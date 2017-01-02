package com.hartwig.hmftools.healthchecker.context;

import com.hartwig.hmftools.common.exception.HartwigException;

import org.jetbrains.annotations.NotNull;

class MalformedRunDirException extends HartwigException {

    private static final String MALFORMED_RUNDIR = "RunDir %s does not match with expectations";

    MalformedRunDirException(@NotNull final String runDirectory) {
        super(String.format(MALFORMED_RUNDIR, runDirectory));
    }
}
