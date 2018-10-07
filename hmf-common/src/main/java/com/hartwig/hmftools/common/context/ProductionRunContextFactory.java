package com.hartwig.hmftools.common.context;

import com.hartwig.hmftools.common.io.exception.MalformedFileException;

import org.jetbrains.annotations.NotNull;

public final class ProductionRunContextFactory {

    private ProductionRunContextFactory() {
    }

    @NotNull
    public static RunContext fromRunDirectory(@NotNull final String runDirectory) throws MalformedFileException {
        final RunContext runContextFromMetaData = MetaDataResolver.fromMetaDataFile(runDirectory);
        if (runContextFromMetaData == null) {
            throw new MalformedFileException("Could not resolve run context from meta data!");
        }
        return runContextFromMetaData;
    }
}
