package com.hartwig.hmftools.patientdb.clinical.context;

import java.io.IOException;

import com.hartwig.hmftools.common.runcontext.MetaDataResolver;
import com.hartwig.hmftools.common.runcontext.RunContext;

import org.jetbrains.annotations.NotNull;

public final class ProductionRunContextFactory {

    private ProductionRunContextFactory() {
    }

    @NotNull
    public static RunContext fromRunDirectory(@NotNull String runDirectory, @NotNull String pipelineVersionFile) throws IOException {
        RunContext runContextFromMetaData = MetaDataResolver.fromMetaDataFile(runDirectory, pipelineVersionFile);
        if (runContextFromMetaData == null) {
            throw new IOException("Could not resolve run context from meta data for " + runDirectory);
        }
        return runContextFromMetaData;
    }
}
