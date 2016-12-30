package com.hartwig.healthchecker.exception;

import org.jetbrains.annotations.NotNull;

public class MalformedRunDirException extends HealthChecksException {

    private static final String MALFORMED_RUNDIR = "RunDir %s does not match with expectations";

    public MalformedRunDirException(@NotNull final String runDirectory) {
        super(String.format(MALFORMED_RUNDIR, runDirectory));
    }
}
