package com.hartwig.hmftools.healthchecker.context;

import com.hartwig.hmftools.common.exception.HartwigException;

import org.jetbrains.annotations.NotNull;

public final class ProductionRunContextFactory {

    private ProductionRunContextFactory() {
    }

    @NotNull
    public static RunContext fromRunDirectory(@NotNull final String runDirectory) throws HartwigException {
        return RunDirectoryResolver.fromRunDirectory(runDirectory);
    }
}
