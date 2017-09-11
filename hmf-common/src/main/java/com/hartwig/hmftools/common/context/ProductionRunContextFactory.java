package com.hartwig.hmftools.common.context;

import com.hartwig.hmftools.common.exception.HartwigException;

import org.jetbrains.annotations.NotNull;

public enum ProductionRunContextFactory {

    ;

    @NotNull
    public static RunContext fromRunDirectory(@NotNull final String runDirectory) throws HartwigException {
        final RunContext runContextFromMetaData = MetaDataResolver.fromMetaDataFile(runDirectory);
        if (runContextFromMetaData == null) {
            throw new HartwigException("Could not resolve run context from meta data!");
        }
        return runContextFromMetaData;
    }
}
