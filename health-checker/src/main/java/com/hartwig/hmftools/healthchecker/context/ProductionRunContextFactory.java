package com.hartwig.hmftools.healthchecker.context;

import com.hartwig.hmftools.common.exception.HartwigException;

import org.jetbrains.annotations.NotNull;

public final class ProductionRunContextFactory {

    private ProductionRunContextFactory() {
    }

    @NotNull
    public static RunContext fromRunDirectory(@NotNull final String runDirectory) throws HartwigException {
        final RunContext runContextFromMetaData = MetaDataResolver.fromMetaDataFile(runDirectory);
        // KODU: Modern runs define their meta data in a specific metadata file, but this is not present for every old
        // run and it would be nice to be able to run health checker on both new and old runs.
        return runContextFromMetaData != null ?
                runContextFromMetaData :
                RunDirectoryResolver.fromRunDirectory(runDirectory);
    }
}
