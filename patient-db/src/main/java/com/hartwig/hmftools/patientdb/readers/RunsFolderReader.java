package com.hartwig.hmftools.patientdb.readers;

import java.io.File;
import java.io.IOException;
import java.text.ParseException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.io.FolderChecker;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class RunsFolderReader {
    private static final Logger LOGGER = LogManager.getLogger(RunsFolderReader.class);

    @NotNull
    public static List<RunContext> getRunContexts(@NotNull final File dir) throws IOException, ParseException {
        final List<RunContext> runContexts = Lists.newArrayList();
        final File[] folders = dir.listFiles(File::isDirectory);
        if (folders == null) {
            throw new IOException("List files in " + dir.getName() + " returned null.");
        }
        for (final File folder : folders) {
            try {
                final String runDirectory = FolderChecker.build().checkFolder(folder.getPath());
                final RunContext runContext = ProductionRunContextFactory.fromRunDirectory(runDirectory);
                runContexts.add(runContext);
            } catch (IOException | HartwigException e) {
                LOGGER.error(e.getMessage());
            }
        }
        return runContexts;
    }
}
