package com.hartwig.hmftools.patientdb.clinical.readers;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.clinical.context.ProductionRunContextFactory;
import com.hartwig.hmftools.patientdb.clinical.context.RunContext;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class RunsFolderReader {

    private static final Logger LOGGER = LogManager.getLogger(RunsFolderReader.class);

    private RunsFolderReader() {
    }

    @NotNull
    public static List<RunContext> extractRunContexts(@NotNull File dir, @NotNull String pipelineVersionFile) throws IOException {
        List<RunContext> runContexts = Lists.newArrayList();
        File[] folders = dir.listFiles(File::isDirectory);
        if (folders == null) {
            throw new IOException("List files in " + dir.getName() + " returned null.");
        }

        for (File folder : folders) {
            if (folder.exists() && folder.isDirectory()) {
                runContexts.add(ProductionRunContextFactory.fromRunDirectory(folder.getPath(), pipelineVersionFile));
            } else {
                LOGGER.warn("Could not process run since file '{}' doesn't seem to be a folder", folder.getPath());
            }
        }
        return runContexts;
    }
}
