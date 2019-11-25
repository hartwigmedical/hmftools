package com.hartwig.hmftools.patientdb.context;

import com.hartwig.hmftools.common.utils.io.exception.MalformedFileException;

import org.jetbrains.annotations.NotNull;

public final class ProductionRunContextFactory {

    private ProductionRunContextFactory() {
    }

    @NotNull
    public static RunContext fromRunDirectory(@NotNull final String runDirectory, @NotNull String whichPackages) throws MalformedFileException {
        final RunContext runContextFromMetaData = MetaDataResolver.fromMetaDataFile(runDirectory, whichPackages);
        if (runContextFromMetaData == null) {
            throw new MalformedFileException("Could not resolve run context from meta data!");
        }
        return runContextFromMetaData;
    }
}
