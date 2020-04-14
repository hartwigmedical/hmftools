package com.hartwig.hmftools.patientdb.context;

import com.hartwig.hmftools.common.utils.io.exception.MalformedFileException;

import org.jetbrains.annotations.NotNull;

public final class ProductionRunContextFactory {

    private ProductionRunContextFactory() {
    }

    @NotNull
    public static RunContext fromRunDirectory(@NotNull String runDirectory) throws MalformedFileException {
        RunContext runContextFromMetaData = MetaDataResolver.fromMetaDataFile(runDirectory);
        if (runContextFromMetaData == null) {
            throw new MalformedFileException("Could not resolve run context from meta data!");
        }
        return runContextFromMetaData;
    }
}
